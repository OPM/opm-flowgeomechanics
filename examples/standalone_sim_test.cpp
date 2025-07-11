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

#include <dune/common/exceptions.hh>
#include <dune/common/filledarray.hh> // needed for printSparseMatrix??
#include <dune/common/indices.hh> // needed for _0, _1, etc.

#include <dune/grid/common/mcmgmapper.hh> // mapper class
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/uggrid.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <dune/istl/io.hh> // needed for printSparseMatrix??
#include <dune/istl/matrixmarket.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/solvers.hh>

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
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace
{
template <class M, class X, class Y>
class ReducedMatrixAdapter; // forward declaration

using Dune::Indices::_0;
using Dune::Indices::_1;

using Grid = Dune::FoamGrid<2, 3>;
using Vector = Dune::BlockVector<double>;
using VectorHP = Dune::MultiTypeBlockVector<Vector, Vector>;
using HTrans = std::tuple<std::size_t, std::size_t, double, double>;

using FullMatrix = Dune::DynamicMatrix<double>;
using SparseMatrix = Dune::BCRSMatrix<double>;
using EquationSystem = std::tuple<std::shared_ptr<FullMatrix>, // aperture matrix
                                  std::shared_ptr<SparseMatrix>, // pressure matrix
                                  std::vector<HTrans>, // inter-cell transmissibility factors
                                  std::vector<double>>; // leakoff factors
using SystemMatrix = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<FullMatrix, SparseMatrix>,
                                                Dune::MultiTypeBlockVector<SparseMatrix, SparseMatrix>>;
using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<Grid::LeafGridView>;
using RMAdapter = ReducedMatrixAdapter<SparseMatrix, Vector, Vector>;

// ----------------------------------------------------------------------------
std::unique_ptr<const Grid>
readgrid(const char* const name)
// ----------------------------------------------------------------------------
{
    try {
        return Dune::GmshReader<Grid>::read(name); // unique_ptr
    } catch (const Dune::Exception& e) {
        std::cerr << "Unable to read file: " << name << std::endl;
    }

    return {};
}

#if 0
// ----------------------------------------------------------------------------
void
dump_matrix(const FullMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    const std::size_t rows = m.N();
    const std::size_t cols = m.M();

    std::ofstream os(name);
    for (std::size_t i = 0; i != rows; ++i) {
        for (std::size_t j = 0; j != cols; ++j) {
            os << m[i][j] << ((j == cols - 1) ? "\n" : " ");
        }
    }
}

// ----------------------------------------------------------------------------
void
dump_matrix(const SparseMatrix& m, const char* const name)
// ----------------------------------------------------------------------------
{
    std::ofstream os(name);

    for (auto rowIt = m.begin(); rowIt != m.end(); ++rowIt) {
        const std::size_t i = rowIt.index();
        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const std::size_t j = colIt.index();
            os << i + 1 << " " << j + 1 << " " << m[i][j] << "\n";
        }
    }
}

// ----------------------------------------------------------------------------
void
dump_vector(const Vector& v, const char* const name)
// ----------------------------------------------------------------------------
{
    std::ofstream os(name);
    for (std::size_t i = 0; i != v.size(); ++i) {
        os << v[i] << "\n";
    }
}
#endif

// ----------------------------------------------------------------------------
template <int N>
std::array<std::size_t, N>
n_closest(const Grid& grid)
// ----------------------------------------------------------------------------
{
    using Elem = std::tuple<std::size_t, double>;
    std::vector<Elem> distances;

    const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());
    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int eIdx = mapper.index(element);
        const auto center = element.geometry().center();
        distances.emplace_back(eIdx, center.two_norm());
    }

    std::sort(distances.begin(), distances.end(), [](Elem& a, Elem& b) {
        return std::get<1>(a) < std::get<1>(b);
    });

    std::array<std::size_t, N> result;
    for (int i = 0; i != N; ++i) {
        result[i] = std::get<0>(distances[i]);
    }

    return result;
}

// ----------------------------------------------------------------------------
std::shared_ptr<FullMatrix>
computeApertureMatrix(const Grid& G, const double young, const double poisson)
// ----------------------------------------------------------------------------
{
    const int nc = G.leafGridView().size(0);

    auto result = std::make_shared<FullMatrix>(nc, nc);
    ddm::assembleMatrix(*result, young, poisson, G);

    *result *= -1;

    return result;
}

// ----------------------------------------------------------------------------
std::tuple<std::vector<HTrans>, std::vector<double>>
computeHTrans(const Grid& grid, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
    std::vector<HTrans> htransvec;
    std::vector<double> leakvec(grid.leafGridView().size(0), 0.0);

    const ElementMapper mapper(grid.leafGridView(), Dune::mcmgElementLayout());

    for (const auto& element : Dune::elements(grid.leafGridView())) {
        const int eIdx = mapper.index(element); // element index
        const auto geom = element.geometry();

        for (const auto& is : Dune::intersections(grid.leafGridView(), element)) {
            if (is.boundary()) {
                continue;
            }

            const int nIdx = mapper.index(is.outside()); // neighbor element index
            const auto igeom = is.geometry();
            if (eIdx < nIdx) {
                const auto iscenter = igeom.center();
                const auto ncenter = is.outside().geometry().center() - iscenter;
                const auto ecenter = geom.center() - iscenter;
                const double area = igeom.volume();
                const double h1 = area / ecenter.two_norm();
                const double h2 = area / ncenter.two_norm();

                htransvec.emplace_back(nIdx, eIdx, h1, h2);
            }
        }

        leakvec[eIdx] = geom.volume() * leakoff_fac;
    }

    return {htransvec, leakvec};
}

// ----------------------------------------------------------------------------
std::tuple<std::shared_ptr<SparseMatrix>, std::vector<HTrans>, std::vector<double>>
computePressureMatrix(const Grid& G, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
    auto mat = std::make_shared<SparseMatrix>();

    // initialize the sparsity pattern of a pressure matrix, with all entries set to zero.
    auto hvec = computeHTrans(G, leakoff_fac);

    const std::size_t nc = G.leafGridView().size(0);

    mat->setBuildMode(SparseMatrix::implicit);
    mat->setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
    mat->setSize(nc, nc);

    for (const auto& elem : std::get<0>(hvec)) {
        const auto ix = static_cast<std::size_t>(std::get<0>(elem));
        const auto jx = static_cast<std::size_t>(std::get<1>(elem));

        mat->entry(ix, jx) = 0.0;
        mat->entry(jx, ix) = 0.0;
        mat->entry(jx, jx) = 0.0;
        mat->entry(ix, ix) = 0.0;
    }

    mat->compress();

    return std::tuple_cat(std::make_tuple(mat), hvec);
}

// ----------------------------------------------------------------------------
EquationSystem
computeEquationSystem(const Grid& G, const double young, const double poisson, const double leakoff_fac)
// ----------------------------------------------------------------------------
{
    return std::tuple_cat(std::make_tuple(computeApertureMatrix(G, young, poisson)),
                          computePressureMatrix(G, leakoff_fac));
}

// ----------------------------------------------------------------------------
template <class M, class X, class Y>
class ReducedMatrixAdapter : public Dune::LinearOperator<X, Y>
// ----------------------------------------------------------------------------
{
public:
    ReducedMatrixAdapter(const M& mat, const std::vector<std::size_t> elim)
        : mat_(mat)
        , elim_(sorted_(elim))
        , keep_(nonelim_(elim_, mat.M()))
        , tmpX_(mat.M())
        , tmpY_(mat.M())
    {
    }

    virtual ~ReducedMatrixAdapter()
    {
    }

    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

    void apply(const X& x, Y& y) const override
    {
        expand(x, tmpX_); // expand x into tmpX_
        mat_.mv(tmpX_, tmpY_);
        contract(tmpY_, y); // contract tmpY_ into y
    }

    void applyscaleadd(field_type alpha, const X& x, Y& y) const override
    {
        expand(x, tmpX_); // expand x into tmpX_
        expand(y, tmpY_); // expand y into tmpY_
        mat_.usmv(alpha, tmpX_, tmpY_);
        contract(tmpY_, y); // contract tmpY_ into y
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

    void expand(const X& source, X& target) const
    {
        target.resize(fullSize());
        target = 0;
        for (std::size_t i = 0; i != reducedSize(); ++i) {
            target[keep_[i]] = source[i];
        }
    }

    void contract(const Y& source, Y& target) const
    {
        target.resize(reducedSize());

        for (std::size_t i = 0; i != reducedSize(); ++i) {
            target[i] = source[keep_[i]];
        }
    }

    std::size_t fullSize() const
    {
        return mat_.M();
    }

    std::size_t reducedSize() const
    {
        return keep_.size();
    }

private:
    std::vector<std::size_t> sorted_(const std::vector<std::size_t>& vec)
    {
        std::vector<std::size_t> result(vec);
        std::sort(result.begin(), result.end());
        return result;
    }

    std::vector<std::size_t> nonelim_(const std::vector<std::size_t>& elim, std::size_t nc)
    {
        // return vector with all incides that are _not_ found in elim
        std::vector<std::size_t> all_ixs(nc);
        std::iota(all_ixs.begin(), all_ixs.end(), 0);

        std::vector<std::size_t> nonelim;
        std::set_difference(
            all_ixs.begin(), all_ixs.end(), elim.begin(), elim.end(), std::back_inserter(nonelim));

        return nonelim;
    }

    const M& mat_;
    const std::vector<std::size_t> elim_;
    const std::vector<std::size_t> keep_;
    mutable X tmpX_;
    mutable Y tmpY_;
};

// ----------------------------------------------------------------------------
template <typename T>
T
hmean(const T a, const T b)
// ----------------------------------------------------------------------------
{
    T tmp = T(1) / (T(1) / a + T(1) / b);
    return std::isnan(tmp) ? 0.0 : tmp;
}

#if 0
// ----------------------------------------------------------------------------
void
updateTrans(SparseMatrix& mat,
            const std::vector<HTrans>& htransvec,
            const Vector& h,
            const std::vector<double>& leakvec)
// ----------------------------------------------------------------------------
{
    mat = 0; // reset matrix before starting to fill in values

    for (const auto& [i, j, t1, t2] : htransvec) {
        assert(i != j);
        const double h1 = h[i];
        const double h2 = h[j];

        const double trans1 = h1 * h1 * t1 / 12.0;
        const double trans2 = h2 * h2 * t2 / 12.0;

        const double T = hmean(trans1, trans2);

        // update diagonal
        mat[i][i] += T;
        mat[j][j] += T;

        // update off-diagonal terms
        mat[i][j] -= T;
        mat[j][i] -= T;
    }

    // add-in leakage term on diagonal
    for (std::size_t ix = 0; ix != mat.N(); ++ix) {
        mat[ix][ix] += leakvec[ix];
    }
}
#endif

// ----------------------------------------------------------------------------
// Note: 'elim' and 'keep' should both be sorted, complementary, and together
// contain all the indices from 0 to M.N() and/or M.M().
// The matrix 'M' will be reduced, and any columns extracted will be collected and
// retuerned. (Rows will just be removed).
SparseMatrix
reduceMatrix(SparseMatrix& M,
             const bool rows,
             const bool cols,
             const std::vector<std::size_t>& elim,
             const std::vector<std::size_t>& keep)
// ----------------------------------------------------------------------------
{
    if (cols) {
        assert(elim.size() + keep.size() == M.M());
    }

    if (rows) {
        assert(elim.size() + keep.size() == M.N());
    }

    SparseMatrix Mreduced, reduced_cols;

    Mreduced.setBuildMode(SparseMatrix::implicit);
    Mreduced.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
    Mreduced.setSize(rows ? keep.size() : M.N(), cols ? keep.size() : M.M());

    reduced_cols.setBuildMode(SparseMatrix::implicit);
    reduced_cols.setImplicitBuildModeParameters(3 * 6 - 3 - 2, 0.4);
    reduced_cols.setSize(Mreduced.N(), cols ? elim.size() : 0);

    std::vector<int> mapping(keep.size() + elim.size(), 0);
    for (std::size_t i = 0; i != keep.size(); ++i) {
        mapping[keep[i]] = i;
    }

    for (std::size_t i = 0; i != elim.size(); ++i) {
        mapping[elim[i]] = -(i + 1);
    }

    for (auto rowIt = M.begin(); rowIt != M.end(); ++rowIt) {
        const std::size_t i = rowIt.index();
        if (rows && mapping[i] < 0) {
            continue; // this row should be eliminated, nothing further to do
        }

        for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
            const std::size_t j = colIt.index();

            if (cols && mapping[j] < 0) {
                // this column should be eliminated
                reduced_cols.entry(rows ? mapping[i] : i, -mapping[i] - 1) = M[i][j];
            } else {
                // this column should be kept
                Mreduced.entry(rows ? mapping[i] : i, cols ? mapping[j] : j) = M[i][j];
            }
        }
    }

    Mreduced.compress();
    reduced_cols.compress();

    std::swap(M, Mreduced);

    return reduced_cols;
}

// ----------------------------------------------------------------------------
class FlowSystemMatrices
// ----------------------------------------------------------------------------
{
public:
    FlowSystemMatrices(const SparseMatrix& M,
                       const std::vector<HTrans>& htransvec,
                       const std::vector<double>& leakvec)
        : M_(M)
        , C_(M)
        , rhs_()
        , htransvec_(htransvec)
        , leakvec_(leakvec)
        , eliminated_indices_()
        , kept_indices_(M.N())
        , remapping_()
        , fixed_bhp_(0)
        , fixed_rate_(0)
        , ratecells_()
    {
        std::iota(kept_indices_.begin(), kept_indices_.end(), 0);
        std::iota(remapping_.begin(), remapping_.end(), 0);
    }

    void setBHPCells(const std::vector<std::size_t>& bhp_cells, const double bhp)
    {
        if (!eliminated_indices_.empty()) {
            throw std::runtime_error("System has already been reduced.  Cannot re-set BHP cells");
        }

        eliminated_indices_ = sorted_(bhp_cells);
        kept_indices_ = noneliminated_(eliminated_indices_, M_.M());
        remapping_ = compute_remapping_(eliminated_indices_, kept_indices_, M_.M());

        reduceMatrix(M_, true, true, eliminated_indices_, kept_indices_);
        reduceMatrix(C_, true, false, eliminated_indices_, kept_indices_);

        rhs_.resize(M_.N());
        fixed_bhp_ = bhp;
    }

    void setRateCells(const std::vector<std::size_t>& rate_cells, const double rate)
    {
        fixed_rate_ = rate;
        ratecells_ = rate_cells;
    }

    std::size_t pnum() const
    {
        return M_.N();
    } // number of (remaining) pressure values

    std::size_t hnum() const
    {
        return C_.N();
    } // number of aperture values

    void update(const VectorHP& hp)
    {
        update_M_and_rhs_(hp);
        update_C_(hp);
    }

    const SparseMatrix& M() const
    {
        return M_;
    }

    const SparseMatrix& C() const
    {
        return C_;
    }

    const Vector& rhs() const
    {
        return rhs_;
    } // should only be called after update

    Vector reduceP(const Vector& p)
    {
        Vector result(pnum());

        result = 0;
        for (std::size_t i = 0; i != pnum(); ++i) {
            result[i] = p[kept_indices_[i]];
        }

        return result;
    }

    Vector expandP(const Vector& p, bool insert_bhp = false)
    {
        assert(p.size() == kept_indices_.size());

        Vector result(remapping_.size());

        result = insert_bhp ? fixed_bhp_ : 0;
        for (std::size_t i = 0; i != p.size(); ++i) {
            result[kept_indices_[i]] = p[i];
        }

        return result;
    }

    const std::vector<std::size_t>& elim()
    {
        return eliminated_indices_;
    }

    const std::vector<std::size_t>& kept()
    {
        return kept_indices_;
    }

    double fixedBHP()
    {
        return fixed_bhp_;
    }

private:
    static std::vector<std::size_t> sorted_(const std::vector<std::size_t>& vec)
    {
        std::vector<std::size_t> result(vec);
        std::sort(result.begin(), result.end());
        return result;
    }

    static std::vector<std::size_t> noneliminated_(const std::vector<std::size_t>& elim,
                                                   const std::size_t nc)
    {
        // return vector with all incides that are _not_ found in elim
        std::vector<std::size_t> all_ixs(nc);
        std::iota(all_ixs.begin(), all_ixs.end(), 0);

        std::vector<std::size_t> nonelim;
        std::set_difference(
            all_ixs.begin(), all_ixs.end(), elim.begin(), elim.end(), std::back_inserter(nonelim));

        return nonelim;
    }

    static std::vector<int> compute_remapping_(const std::vector<std::size_t>& eliminated,
                                               const std::vector<std::size_t>& kept,
                                               const std::size_t nc)
    {
        std::vector<int> result(nc);
        for (std::size_t i = 0; i != eliminated.size(); ++i) {
            result[eliminated[i]] = -i - 1;
        }

        for (std::size_t i = 0; i != kept.size(); ++i) {
            result[kept[i]] = i;
        }

        return result;
    }

    void update_M_and_rhs_(const VectorHP& hp)
    {
        // update M_ and rhs_reduct_
        M_ = 0;
        rhs_ = 0;

        for (const auto& [i, j, t1, t2] : htransvec_) {
            assert(i != j);

            const double h1 = hp[_0][i];
            const double h2 = hp[_0][j];

            const double trans1 = h1 * h1 * t1 / 12.0;
            const double trans2 = h2 * h2 * t2 / 12.0;

            const double T = hmean(trans1, trans2);

            // remapped indices
            const int ir = remapping_[i];
            const int jr = remapping_[j];

            // update diagonal
            if (ir >= 0) {
                M_[ir][ir] += T;
            }

            if (jr >= 0) {
                M_[jr][jr] += T;
            }

            // update off-diagonal terms, and the reduction terms of rhs_reduct_
            if (ir >= 0) {
                if (jr >= 0) {
                    M_[ir][jr] -= T;
                } else {
                    rhs_[ir] += T * fixed_bhp_;
                }
            }

            if (jr >= 0) {
                if (ir >= 0) {
                    M_[jr][ir] -= T;
                } else {
                    rhs_[jr] += T * fixed_bhp_;
                }
            }
        }

        // add leakage terms
        for (std::size_t ix = 0; ix != leakvec_.size(); ++ix) {
            if (remapping_[ix] >= 0) {
                M_[remapping_[ix]][remapping_[ix]] += leakvec_[ix];
            }
        }

        // add source terms
        for (const auto& ix : ratecells_) {
            if (remapping_[ix] >= 0) {
                rhs_[remapping_[ix]] = fixed_rate_;
            }
        }
    }

    void update_C_(const VectorHP& hp)
    {
        C_ = 0;
        // intermediary matrices for computing the partial derivatives
        for (const auto& [i, j, t1, t2] : htransvec_) {
            assert(i != j);

            const int ir = remapping_[i];
            const int jr = remapping_[j];

            const double h1 = hp[_0][i];
            const double h2 = hp[_0][j];
            const double p1 = hp[_1][i];
            const double p2 = hp[_1][j];

            const double q = h1 * h1 * h2 * h2 * t1 * t2; // numerator
            const double d1q = 2 * h1 * h2 * h2 * t1 * t2;
            const double d2q = h1 * h1 * 2 * h2 * t1 * t2;

            const double r = 12 * (h1 * h1 * t1 + h2 * h2 * t2); // denominator
            const double d1r = 24 * h1 * t1;
            const double d2r = 24 * h2 * t2;

            const double dTdh1 = (r == 0) ? 0.0 : (d1q * r - q * d1r) / (r * r);
            const double dTdh2 = (r == 0) ? 0.0 : (d2q * r - q * d2r) / (r * r);

            // diagonal elements
            if (ir >= 0) {
                C_[ir][i] += dTdh1 * (p1 - p2);
            }

            if (jr >= 0) {
                C_[jr][j] += dTdh2 * (p2 - p1);
            }

            // off-diagonal elements
            if (ir >= 0) {
                C_[ir][j] += dTdh2 * (p1 - p2);
            }

            if (jr >= 0) {
                C_[jr][i] += dTdh1 * (p2 - p1);
            }
        }
    }

    SparseMatrix M_;
    SparseMatrix C_;
    Vector rhs_;
    const std::vector<HTrans>& htransvec_;
    const std::vector<double>& leakvec_;
    std::vector<std::size_t> eliminated_indices_;
    std::vector<std::size_t> kept_indices_;
    std::vector<int> remapping_;
    double fixed_bhp_;
    double fixed_rate_;
    std::vector<std::size_t> ratecells_;
};

// ----------------------------------------------------------------------------
SparseMatrix
makeIdentity(std::size_t num, double fac = 1)
// ----------------------------------------------------------------------------
{
    // build a sparse matrix to represent the identity matrix of a given size
    SparseMatrix M;

    M.setBuildMode(SparseMatrix::implicit);
    M.setImplicitBuildModeParameters(1, 0.1);
    M.setSize(num, num);

    for (std::size_t i = 0; i != num; ++i) {
        M.entry(i, i) = fac;
    }

    M.compress();

    return M;
}

// ----------------------------------------------------------------------------
template <typename Mat>
Vector
diagvec(const Mat& M)
// ----------------------------------------------------------------------------
{
    Vector res(M.M());

    for (std::size_t i = 0; i != res.size(); ++i) {
        res[i] = M[i][i];
    }

    return res;
}

// ----------------------------------------------------------------------------
class TailoredPrecondFull : public Dune::Preconditioner<VectorHP, VectorHP>
// ----------------------------------------------------------------------------
{
public:
    TailoredPrecondFull(const SystemMatrix& S)
        : A_(S[_0][_0])
        , M_diag_(diagvec(S[_1][_1]))
    {
    }

    void apply(VectorHP& v, const VectorHP& d) override
    {
        A_.solve(v[_0], d[_0]);

        for (std::size_t i = 0; i != M_diag_.size(); ++i) {
            v[_1][i] = d[_1][i] / M_diag_[i];
        }
    }

    void post([[maybe_unused]] VectorHP& v) override
    {
    }

    void pre([[maybe_unused]] VectorHP& x, [[maybe_unused]] VectorHP& b) override
    {
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const FullMatrix& A_;
    const Vector M_diag_;
};

// ----------------------------------------------------------------------------
class TailoredPrecondDiag : public Dune::Preconditioner<VectorHP, VectorHP>
// ----------------------------------------------------------------------------
{
public:
    TailoredPrecondDiag(const SystemMatrix& S)
        : A_diag_(diagvec(S[_0][_0]))
        , M_diag_(diagvec(S[_1][_1]))
    {
    }

    void apply(VectorHP& v, const VectorHP& d) override
    {
        for (std::size_t i = 0; i != A_diag_.size(); ++i) {
            v[_0][i] = d[_0][i] / A_diag_[i];
        }

        for (std::size_t i = 0; i != M_diag_.size(); ++i) {
            v[_1][i] = d[_1][i] / M_diag_[i];
        }
    }

    void post([[maybe_unused]] VectorHP& v) override
    {
    }

    void pre([[maybe_unused]] VectorHP& x, [[maybe_unused]] VectorHP& b) override
    {
    }

    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }

private:
    const Vector A_diag_;
    const Vector M_diag_;
};

// ----------------------------------------------------------------------------
double
estimate_step_fac(const VectorHP& x, const VectorHP& dx)
// ----------------------------------------------------------------------------
{
    // estimate what might be a safe step size to avoid exiting the convergence radius

    const double f1 = dx[_0].infinity_norm() / x[_0].infinity_norm();
    const double f2 = dx[_1].infinity_norm() / x[_1].infinity_norm();
    const double fmax = std::max(f1, f2);
    const double threshold = 0.95;
    const double fac_min = 1e-2; // 1e-2;
    const double fac = (fmax < threshold) ? 1.0 : std::max(threshold / fmax, fac_min);

    return fac;
}

// ----------------------------------------------------------------------------
bool
nonlinearIteration(const FullMatrix& A,
                   FlowSystemMatrices& SFM, // will be updated
                   const VectorHP& hp,
                   VectorHP& dhp, // increment (to be computed)
                   const std::function<bool(const VectorHP&)>& converged_pred)
// ----------------------------------------------------------------------------
{
    SFM.update(hp);

    auto negI = makeIdentity(A.N(), -1);
    auto negI_rhs = reduceMatrix(negI, false, true, SFM.elim(), SFM.kept());

    SystemMatrix S {{A, negI}, {SFM.C(), SFM.M()}};
    VectorHP rhs;

    Vector tmp(SFM.elim().size());
    for (std::size_t i = 0; i != SFM.elim().size(); ++i) {
        tmp[i] = SFM.fixedBHP() * -1; // sign flip, since we are moving to right-hand side
    }

    rhs[_0].resize(A.N());
    negI_rhs.mv(tmp, rhs[_0]); // accounting for eliminated pressure values

    rhs[_1] = SFM.rhs();

    dhp = rhs; // ensure it has the right number of elements

    VectorHP hp_reduced = hp;
    hp_reduced[_1] = SFM.reduceP(hp[_1]);

    auto S0 = S;
    S0[_1][_0] = 0; // The C-matrix should be zero here
    S0.mmv(hp_reduced, rhs); // rhs = rhs - S0 * hp_reduced

    std::cout << "Residual norm is: " << rhs.infinity_norm() << std::endl;
    if (converged_pred(rhs)) {
        return true; // system converged
    }

    // solve system equations
    Dune::MatrixAdapter<SystemMatrix, VectorHP, VectorHP> S_linop(S);

    TailoredPrecondDiag precond(S);
    Dune::InverseOperatorResult iores;

    auto psolver = Dune::BiCGSTABSolver<VectorHP>(S_linop,
                                                  precond,
                                                  1e-20, // 1e-15, // desired rhs reduction factor
                                                  200, // max number of iterations
                                                  1); // verbose

    VectorHP res_copy = rhs; // we need to keep the original rhs unmodified
    psolver.apply(dhp, res_copy, iores);

    dhp[_1] = SFM.expandP(dhp[_1], false);

    std::cout << "hp:  " << hp[_0].infinity_norm() << " " << hp[_1].infinity_norm() << std::endl;
    std::cout << "dhp: " << dhp[_0].infinity_norm() << " " << dhp[_1].infinity_norm() << std::endl;

    // the following is a heuristic way to try to limit stepsize to stay within
    // convergence radius
    const double step_fac = estimate_step_fac(hp, dhp);
    std::cout << "fac: " << step_fac << std::endl;
    dhp *= step_fac;

    return false; // system has not yet been shown to have converged
}

void
cap_vals(VectorHP& v, const double low_h, const double low_p, const double high_h, const double high_p)
{
    for (std::size_t i = 0; i != v[_0].size(); ++i) {
        v[_0][i] = std::clamp(v[_0][i], low_h, high_h);
        v[_1][i] = std::clamp(v[_1][i], low_p, high_p);
    }
}

// ----------------------------------------------------------------------------
int
solveCoupledFull(Vector& p,
                 Vector& h,
                 const EquationSystem& eqsys,
                 const std::vector<std::size_t> bhpcells,
                 const double bhp,
                 const std::vector<std::size_t> ratecells,
                 const double rate,
                 const double convergence_tol = 1e-5,
                 const int max_nonlin_iter = 40)
// ----------------------------------------------------------------------------
{
    // We consider the nonlinear system:
    // | A          -I |  | h |   | 0 |   (mechanics equation)
    // |               |  |   | = |   |
    // | C(h, p)  M(h) |  | p |   | q |   ( flow equation)

    const FullMatrix& A = *std::get<0>(eqsys);
    const SparseMatrix& M = *std::get<1>(eqsys);
    const std::vector<HTrans>& htransvec = std::get<2>(eqsys);
    const std::vector<double>& leakvec = std::get<3>(eqsys);

    VectorHP hp {h, p}; // solution vector, grouping aperture 'h' and pressure 'p'

    for (std::size_t i = 0; i != bhpcells.size(); ++i)
        hp[_1][bhpcells[i]] = bhp; // these are the fixed values

    FlowSystemMatrices FSM(M, htransvec, leakvec);
    FSM.setBHPCells(bhpcells, bhp);
    FSM.setRateCells(ratecells, rate);

    std::function<bool(const VectorHP&)> converged_pred = [&](const VectorHP& res) {
        std::cout << res[_0].infinity_norm() << ", " << res[_1].infinity_norm() << std::endl;
        const double aperture_tol = (bhpcells.size() > 0) ? convergence_tol * bhp : convergence_tol;
        const double flow_tol = convergence_tol;
        return res[_0].infinity_norm() < aperture_tol && res[_1].infinity_norm() < flow_tol;
    };

    VectorHP dhp;
    int iter = 0;

    // nonlinear solve loop
    while (!nonlinearIteration(A, FSM, hp, dhp, converged_pred) && ++iter < max_nonlin_iter) {
        hp += dhp;
        cap_vals(hp, 1e-12, 0, 1e99, 1e99);
        std::cout << "Completed iteration: " << iter << std::endl;
    }

    if (iter == max_nonlin_iter) {
        std::cout << "System did not converge in max allowed number of iterations." << std::endl;
        return -1;
    } else {
        std::cout << "Converged in " << iter << " iterations." << std::endl;
    }

    // system converged, unpacking results
    h = hp[_0];
    p = hp[_1];
    return 0;
}

#if 0
// ----------------------------------------------------------------------------
void
solveCoupledSplit(Vector& p,
                  Vector& h,
                  const EquationSystem& eqsys,
                  const std::vector<std::size_t> bhpcells,
                  const double bhp,
                  const std::vector<std::size_t> ratecells,
                  const double rate,
                  const double convergence_tol = 1e-4,
                  const int max_nonlin_iter = 400)
// ----------------------------------------------------------------------------
{
    // NB: there should be no overlap between the cells in 'bhpcells' and 'ratecells'
    // reduce pressure system

    // Dune::MatrixAdapter<FullMatrix, Vector, Vector> MA_h(*std::get<0>(eqsys));
    const FullMatrix& M_h = *std::get<0>(eqsys);
    SparseMatrix& M_p = *std::get<1>(eqsys);
    const std::vector<HTrans>& htransvec = std::get<2>(eqsys);
    const std::vector<double>& leakvec = std::get<3>(eqsys);

    RMAdapter MA_p(M_p,
                   bhpcells); // allows reducing the system without creating new matrices

    // define convergence criterion
    auto converged = [&](const Vector& v1, const Vector& v2) {
        auto diff = v1;
        diff -= v2; // diff = v1 - v2
        const double delta = diff.infinity_norm();
        const double denom = std::max(v1.infinity_norm(), v2.infinity_norm());
        std::cout << "Delta: " << delta << " Denom: " << denom;
        std::cout << " Delta/denom: " << delta / denom << std::endl;
        return delta / denom < convergence_tol;
    };

    // helper function to set vector values
    auto setvals = [](Vector& v, const std::vector<std::size_t>& ixs, const double val) {
        for (const auto& i : ixs) {
            v[i] = val;
        }
    };

    // initialize pressure
    p = 0;
    if (bhpcells.size() > 0) {
        setvals(p, bhpcells, bhp);
    } else {
        p = 1e6; // we need something to get started without h being 0.
    }

    // solve for aperture
    M_h.solve(h, p); // solve for aperture (h) given pressure (p)

    Vector htmp(h.N());
    htmp = 0;
    Vector rhs_full_rate(p.N()), rhs_full(p.N());
    Vector rhs;
    Dune::InverseOperatorResult iores;
    rhs_full_rate = 0;
    setvals(rhs_full_rate, ratecells, rate);

    int i;
    for (i = 0; i != max_nonlin_iter; ++i) {
        // update pressure matrix entries
        updateTrans(M_p, htransvec, h, leakvec);

        // determine right-hand side modifications from eliminated degrees of freedom

        M_p.mv(p, rhs_full); // only imposed values of p should be nonzero here.
        rhs_full *= -1;
        rhs_full += rhs_full_rate;

        MA_p.contract(rhs_full,
                      rhs); // remove entries corresponding to eliminated equations

        // solve for pressure
        Dune::Richardson<Vector, Vector> precond(1); // "no" preconditioner

        auto psolver = Dune::CGSolver<Vector>(MA_p,
                                              precond,
                                              1e-7, // desired residual reduction factor
                                              100, // max number of iterations
                                              1); // verbose

        Vector p_reduced(MA_p.reducedSize());
        psolver.apply(p_reduced, rhs, iores);
        MA_p.expand(p_reduced, p);
        setvals(p, bhpcells, bhp);

        // solve for aperture again
        M_h.solve(htmp, p); // solve for aperture (h) given pressure (p)

        if (converged(h, htmp)) {
            break;
        }

        // h = htmp;
        h *= 0.0;
        htmp *= 1.0;
        h += htmp;
    }

    if (i == max_nonlin_iter) {
        std::cout << "Warning, did not converge in max number of nonlinear iterations." << std::endl;
    } else {
        std::cout << "Converged in: " << i << " iterations." << std::endl;
    }
}
#endif

} // end anonymous namespace

// ============================================================================
int
main(int varnum, char** vararg)
// ============================================================================
{
    // check argument list
    if (varnum < 2) {
        std::cout << "Grid file name must be given as argument." << std::endl;
        return -1;
    }

    // read grid
    const auto grid = readgrid(vararg[1]);
    if (!grid) {
        return -1;
    }

    // identify well cells
    const std::vector<std::size_t> wellcells {0, 1, 2, 3, 4, 5};

    std::cout << "Identified wellcells: ";
    std::copy(wellcells.begin(), wellcells.end(), std::ostream_iterator<std::size_t>(std::cout, " "));
    std::cout << std::endl;

    // compute equation system
    const double young = 1e9; // fYoung's modulus
    const double poisson = 0.25; // Poisson's ratio
    // const double leakoff_fac = 1e-13;// this breaks for rate
    const double leakoff_fac = 1e-7; // 1e-9; //1e-6; //1e-13; // a bit heuristic;
                                     // conceptually rock perm divided by distance
    const auto eqsys = computeEquationSystem(*grid, young, poisson, leakoff_fac);

    // solve fixed pressure system
    const std::size_t nc = grid->leafGridView().size(0);
    Vector pressure(nc), aperture(nc);
    pressure = 0;
    aperture = 1e-2;
    const double bhp = 1e6; // in Pascal
    const double rate = 0.1 / wellcells.size(); // positive value is _injection_

    // solve fixed rate system
    const auto pressure_cells = std::vector<std::size_t>();
    const auto rate_cells = std::vector<std::size_t>(wellcells.begin(), wellcells.end());

    solveCoupledFull(pressure, aperture, eqsys, pressure_cells, bhp, rate_cells, rate, 1e-5, 1000);

    // write output
    Dune::VTKWriter<Grid::LeafGridView> vtkwriter(grid->leafGridView(), Dune::VTK::nonconforming);
    vtkwriter.addCellData(aperture, "aperture");
    vtkwriter.addCellData(pressure, "pressure");
    vtkwriter.write("output");
} // end main
