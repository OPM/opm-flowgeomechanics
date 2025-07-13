/*
  Copyright 2025 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <config.h>

#include <dune/common/filledarray.hh> // needed for printSparseMatrix??

#include <dune/foamgrid/foamgrid.hh>

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/istl/io.hh> // needed for printSparseMatrix??
#include <dune/istl/matrixmarket.hh>

//@@ there must be a more correct way to ensure UMFpack is included here
#define HAVE_SUITESPARSE_UMFPACK 1
#include <dune/istl/umfpack.hh>

#include <opm/simulators/linalg/FlexibleSolver.hpp>
#include <opm/simulators/linalg/FlowLinearSolverParameters.hpp>
#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/setupPropertyTree.hpp>

#include <opm/geomech/DiscreteDisplacement.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace
{
// definitions and constants
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;
using Grid = Dune::FoamGrid<2, 3>;
using Htrans = std::tuple<std::size_t, std::size_t, double, double>;
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
using PressureOperatorType = Dune::MatrixAdapter<Matrix, Vector, Vector>;
using FlexibleSolverType = Dune::FlexibleSolver<PressureOperatorType>;
using PressureMatrixInfo = std::tuple<std::unique_ptr<Matrix>, std::vector<Htrans>>;

using IntFloatPair = std::tuple<int, double>;

constexpr auto E = 1e9;
constexpr auto nu = 0.25;
constexpr auto init_frac_press = 1e6;

// ============================================================================
std::vector<Htrans>
computeHtrans(const Grid& grid)
// ============================================================================
{
    std::vector<Htrans> result {};

    const auto mapper = ElementMapper {grid.leafGridView(), Dune::mcmgElementLayout()};
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto eIdx = mapper.index(element);
        const auto geom = element.geometry();

        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.boundary()) {
                continue;
            }

            const auto nIdx = mapper.index(is.outside());
            if (nIdx <= eIdx) {
                continue;
            }

            const auto& igeom = is.geometry();

            const auto iscenter = igeom.center();
            const auto ncenter = is.outside().geometry().center() - iscenter;
            const auto ecenter = geom.center() - iscenter;

            const double area = igeom.volume();
            const double h1 = ecenter.two_norm() / area;
            const double h2 = ncenter.two_norm() / area;

            result.emplace_back(nIdx, eIdx, h1, h2);
        }
    }

    return result;
}

// ============================================================================
PressureMatrixInfo
pressureMatrixStructure(const Grid& grid)
// ============================================================================
{
    // computing pre-transmissibilities
    const std::vector<Htrans> htrans = computeHtrans(grid);

    // setting up matrix and initializing its sparsity structure based on the
    // computed pre-transmissibilities
    auto pressure_matrix = std::make_unique<Matrix>();
    pressure_matrix->setBuildMode(Matrix::implicit);

    // map from dof=3*nodes at a cell (ca 3*3*3) to cell
    const auto nc = grid.leafGridView().size(0);

    pressure_matrix->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
    pressure_matrix->setSize(nc, nc);

    for (const auto& elem : htrans) {
        const std::size_t ix = std::get<0>(elem);
        const std::size_t jx = std::get<1>(elem);

        pressure_matrix->entry(ix, jx) = 0.0;
        pressure_matrix->entry(jx, ix) = 0.0;
        pressure_matrix->entry(jx, jx) = 0.0;
        pressure_matrix->entry(ix, ix) = 0.0;
    }

    pressure_matrix->compress();

    return PressureMatrixInfo {std::move(pressure_matrix), htrans};
}

// ============================================================================
template <typename T>
T
hmean(const T a, const T b)
{
    return T(1) / (T(1) / a + T(1) / b);
}
// ============================================================================

// ============================================================================
template <typename T, typename Container>
bool
inrange(T el, Container con)
{
    return std::find(con.begin(), con.end(), el) != con.end();
}
// ============================================================================

// ============================================================================
void
updateTrans(PressureMatrixInfo& pmat,
            const Vector& aperture,
            const std::vector<std::size_t>& imposed_vals_ixs = std::vector<std::size_t>())
// ============================================================================
{
    auto& matrix = *std::get<0>(pmat);
    matrix = 0;

    for (const auto& [i, j, t1, t2] : std::get<1>(pmat)) {
        assert(i != j);

        const double h1 = aperture[i];
        const double h2 = aperture[j];

        const double trans1 = h1 * h1 * t1 / 12.0;
        const double trans2 = h2 * h2 * t2 / 12.0;

        const double T = hmean(trans1, trans2);

        // imposed values will be associated with trivial equations
        const bool i_imposed = inrange(i, imposed_vals_ixs);
        const bool j_imposed = inrange(j, imposed_vals_ixs);

        matrix[i][i] = i_imposed ? 1.0 : static_cast<double>(matrix[i][i]) + T;
        matrix[i][j] = i_imposed ? 0.0 : static_cast<double>(matrix[i][j]) - T;

        matrix[j][j] = j_imposed ? 1.0 : static_cast<double>(matrix[j][j]) + T;
        matrix[j][i] = j_imposed ? 0.0 : static_cast<double>(matrix[j][i]) - T;
    }
}

// ============================================================================
Vector
solvePressure(const Vector& aperture,
              PressureMatrixInfo& pmat,
              const std::vector<IntFloatPair>& fixed_pvals = std::vector<IntFloatPair>())
// ============================================================================
{
    // fill in pressure matrix with actual values, based on current aperture and
    // imposed pressure
    std::vector<std::size_t> imposed_ixs {};
    imposed_ixs.reserve(fixed_pvals.size());

    std::transform(fixed_pvals.begin(),
                   fixed_pvals.end(),
                   std::back_inserter(imposed_ixs),
                   [&](const auto& el) { return std::get<0>(el); });

    updateTrans(pmat, aperture, imposed_ixs);

    // prepare right-hand side
    Vector rhs(aperture.size());
    rhs = 0;
    for (const auto& fp : fixed_pvals) {
        rhs[std::get<0>(fp)] = std::get<1>(fp);
    }

    // setup solver
    const auto prm = Opm::setupPropertyTree(Opm::FlowLinearSolverParameters(), true, true);
    auto op = PressureOperatorType(*std::get<0>(pmat));
    auto psolver_dummy = FlexibleSolverType(op, prm, std::function<Vector()>(), std::size_t {0});

    Dune::UMFPack psolver(*std::get<0>(pmat));

    // solve pressure system
    Vector result(aperture.size());
    Dune::InverseOperatorResult r;
    psolver.apply(result, rhs, r);

    return result;
}

// ============================================================================
template <int N>
std::array<std::size_t, N>
n_closest(const Grid& grid)
// ============================================================================
{
    using Elem = std::tuple<std::size_t, double>;

    std::vector<Elem> distances {};

    const auto mapper = ElementMapper {grid.leafGridView(), Dune::mcmgElementLayout()};

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const auto center = element.geometry().center();
        distances.emplace_back(distances.size(), center.two_norm());
    }

    std::sort(distances.begin(), distances.end(), [](const Elem& a, const Elem& b) {
        return std::get<1>(a) < std::get<1>(b);
    });

    std::array<std::size_t, N> result;
    for (int i = 0; i != N; ++i) {
        result[i] = std::get<0>(distances[i]);
    }

    return result;
}

} // end anonymous namespace

// ============================================================================
int
main()
// ============================================================================
{
    // read test grid from disk
    const auto grid = Dune::GmshReader<Grid>::read("disk.msh"); // unique_ptr

    // compute mech matrix
    const int nc = grid->leafGridView().size(0);

    auto frac_matrix = std::make_unique<Dune::DynamicMatrix<double>>();
    frac_matrix->resize(nc, nc);

    ddm::assembleMatrix(*frac_matrix, E, nu, *grid);
    *frac_matrix *= -1;

    // computing aperture
    Vector frac_press(nc), frac_aperture(nc);
    frac_press = init_frac_press;
    frac_aperture = 0;

    frac_matrix->solve(frac_aperture, frac_press);

    // determining well cells and prescribing fixed pressure at these
    const auto wellcells = n_closest<2>(*grid);

    const double pfixedval = 2.0;

    std::vector<IntFloatPair> fixed_pvals {};
    fixed_pvals.reserve(wellcells.size());
    std::transform(wellcells.begin(),
                   wellcells.end(),
                   std::back_inserter(fixed_pvals),
                   [pfixedval](const std::size_t ix) { return IntFloatPair(ix, pfixedval); });

    // prepare pressure matrix structure
    auto pmat = pressureMatrixStructure(*grid);

    frac_press = solvePressure(frac_aperture, pmat, fixed_pvals);

    std::cout << frac_press;
}
