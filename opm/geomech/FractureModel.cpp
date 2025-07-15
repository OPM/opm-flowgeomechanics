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

#include <opm/geomech/FractureModel.hpp>

#include <dune/common/fmatrixev.hh>

#include <opm/input/eclipse/Schedule/ScheduleState.hpp>
#include <opm/input/eclipse/Schedule/Well/Connection.hpp>
#include <opm/input/eclipse/Schedule/Well/WellConnections.hpp>
#include <opm/input/eclipse/Schedule/Well/WellFractureSeeds.hpp>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/wells/RuntimePerforation.hpp>
#include <opm/simulators/wells/SingleWellState.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <opm/geomech/DiscreteDisplacement.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm
{
PropertyTree
makeDefaultFractureParam()
{
    using namespace std::string_literals;

    PropertyTree fracture_param;
    fracture_param.put("hasfractures", false);
    fracture_param.put("add_perfs_to_schedule", true);
    // solution method
    fracture_param.put("solver.method", "PostSolve"s);
    fracture_param.put("solver.implicit_flow", false);
    fracture_param.put("solver.max_mech_it", 2);

    // solve method for fracture
    fracture_param.put("fractureparam.method.iterate", true);
    fracture_param.put("fractureparam.method.max_it", 3);
    fracture_param.put("fractureparam.method.tolerance", 1e-4);

    //
    fracture_param.put("fractureparam.reduce_boundary", false);
    fracture_param.put("fractureparam.addconnections", true);

    // very experimental to calculate stress contributions from fracture to cell values
    fracture_param.put("fractureparam.include_fracture_contributions", false);

    // seed to be in input file
    fracture_param.put("fractureparam.config.type", "well_seed"s);
    fracture_param.put("fractureparam.config.initial_fracture_width", 1e-4);
    fracture_param.put("fractureparam.config.min_width", 0.0);
    fracture_param.put("fractureparam.config.trires", 10);

    fracture_param.put("fractureparam.solver.method", "if_propagate_trimesh"s);
    fracture_param.put("fractureparam.solver.target_cellcount", 100);
    fracture_param.put("fractureparam.solver.cellcount_threshold", 400);
    fracture_param.put("fractureparam.solver.numcell_threshold", 50);
    fracture_param.put("fractureparam.solver.max_num_coarsening", 200);
    fracture_param.put("fractureparam.solver.efac", 0.5);
    fracture_param.put("fractureparam.solver.rfac", 0.1);
    fracture_param.put("fractureparam.solver.max_expand_iter", 20);
    fracture_param.put("fractureparam.solver.max_iter", 100);
    fracture_param.put("fractureparam.solver.tolerance", 1e-6);
    fracture_param.put("fractureparam.solver.damping", 1e0);
    fracture_param.put("fractureparam.solver.min_width", 1e-3);
    fracture_param.put("fractureparam.solver.max_width", 0.5);
    fracture_param.put("fractureparam.solver.max_dwidth", 5e-3);
    fracture_param.put("fractureparam.solver.max_dp", 1e9);
    fracture_param.put("fractureparam.solver.max_change", 1e5);
    fracture_param.put("fractureparam.solver.verbosity", 0);

    // fracture linear solve
    fracture_param.put("fractureparam.solver.linsolver.tol", 1e-10);
    fracture_param.put("fractureparam.solver.linsolver.max_iter", 1000);
    fracture_param.put("fractureparam.solver.linsolver.verbosity", 0);

    // reservoir fracture coupling
    fracture_param.put("fractureparam.reservoir.dist", 1e0);
    fracture_param.put("fractureparam.reservoir.calculate_dist", false);
    fracture_param.put("fractureparam.reservoir.mobility", 1.3e-3);
    fracture_param.put("fractureparam.reservoir.perm", 1e-13);

    // well fracture coupling
    fracture_param.put("fractureparam.control.type", "perf_pressure"s);
    fracture_param.put("fractureparam.control.rate", 2.9e-2);
    fracture_param.put("fractureparam.control.WI", 1.0e-11);

    fracture_param.put("fractureparam.KMax", 1e6); // in input file
    fracture_param.put("fractureparam.extended_fractures", true);
    fracture_param.put("fractureparam.fractureWI", 0.1);
    fracture_param.put("fractureparam.write_pressure_system", false);
    fracture_param.put("fractureparam.write_fracture_system", false);
    fracture_param.put("fractureparam.pressuresolver", "umfpack"s);
    fracture_param.put("fractureparam.fracturesolver", "notused"s);

    return fracture_param;
}

void
FractureModel::addWell(const std::string& name,
                       const std::vector<FractureWell::Connection>& conns,
                       const std::vector<Point3D>& points,
                       const std::vector<Segment>& segments)
{
    const std::string outputdir = prm_.get<std::string>("outputdir");
    const std::string casename = prm_.get<std::string>("casename");

    // add with no fractures
    wells_.emplace_back(outputdir, casename, name, conns, points, segments);
    well_fractures_.emplace_back();
}

void
FractureModel::addFractures(const ScheduleState& sched)
{
    const auto fracture_type = this->prm_.get<std::string>("config.type", "well_seed");

    if (fracture_type == "perp_well") {
        this->addFracturesPerpWell();
    } else if (fracture_type == "tensile_fracture") {
        this->addFracturesTensile();
    } else if (fracture_type == "well_seed") {
        this->addFracturesWellSeed(sched);
    } else {
        OPM_THROW(std::runtime_error, "Fracture type '" + fracture_type + "' is not supported");
    }

    std::cout << "Added fractures to " << wells_.size() << " wells\n"
              << "Total number of fractures_wells: " << well_fractures_.size() << std::endl;

    int count_frac = 0;
    for (auto i = 0 * well_fractures_.size(); i < well_fractures_.size(); ++i) {
        const auto nfrac = well_fractures_[i].size();
        const auto* pl = (nfrac != 1) ? "s" : "";

        count_frac += static_cast<int>(nfrac);

        std::cout << "Well " << wells_[i].name() << " has " << nfrac << " fracture" << pl << '\n';
    }

    std::cout << "Total number of fractures: " << count_frac << std::endl;
}

void
FractureModel::initFractureStates()
{
    for (auto& fractures : this->well_fractures_) {
        for (auto& fracture : fractures) {
            fracture.initFractureStates();
        }
    }
}

// probably this should be collected in one loop
Dune::FieldVector<double, 6>
FractureModel::stress(const Dune::FieldVector<double, 3>& obs) const
{
    Dune::FieldVector<double, 6> stress {};

    for (const auto& fractures : this->well_fractures_) {
        for (const auto& fracture : fractures) {
            stress += fracture.stress(obs);
        }
    }

    return stress;
}

Dune::FieldVector<double, 6>
FractureModel::strain(const Dune::FieldVector<double, 3>& obs) const
{
    Dune::FieldVector<double, 6> strain {};

    for (const auto& fractures : this->well_fractures_) {
        for (const auto& fracture : fractures) {
            strain += fracture.strain(obs);
        }
    }

    return strain;
}

Dune::FieldVector<double, 3>
FractureModel::disp(const Dune::FieldVector<double, 3>& obs) const
{
    Dune::FieldVector<double, 3> disp {};

    for (const auto& fractures : this->well_fractures_) {
        for (const auto& fracture : fractures) {
            disp += fracture.disp(obs);
        }
    }

    return disp;
}

void
FractureModel::write(const int reportStep) const
{
    auto i = std::size_t {0};

    for (const auto& well : this->wells_) {
        well.write();

        for (const auto& fracture : this->well_fractures_[i]) {
            fracture.write(reportStep);
        }

        ++i;
    }
}

void
FractureModel::writemulti(const double time) const
{
    auto i = std::size_t {0};

    for (const auto& well : this->wells_) {
        if (vtkwritewells_) {
            well.writemulti(time);
        }

        for (const auto& fracture : this->well_fractures_[i]) {
            fracture.writemulti(time);
        }

        ++i;
    }
}

void
FractureModel::updateReservoirProperties()
{
    for (auto& fractures : this->well_fractures_) {
        for (auto& fracture : fractures) {
            fracture.updateReservoirProperties();
        }
    }
}

std::vector<RuntimePerforation>
FractureModel::getExtraWellIndices(const std::string& wellname) const
{
    if (!prm_.get<bool>("addconnections")) {
        return {};
    }

    // for now just do a search
    int well_idx = -1;
    for (std::size_t i = 0; i < wells_.size(); ++i) {
        if (!wells_[i].isActive()) {
            if (wells_[i].name() == wellname) {
                well_idx = i;
            }
            continue;
        }

        if (wells_[i].name() == wellname) {
            well_idx = i;
            // collect all from a well
            std::vector<RuntimePerforation> wellindices;
            for (const auto& frac : well_fractures_[i]) {
                if (!frac.isActive()) {
                    continue;
                }
                const auto perfs = frac.wellIndices();
                wellindices.insert(wellindices.end(), perfs.begin(), perfs.end());
            }

            return wellindices;
        }
    }

    if (false) {
        std::cout << "Well " << wellname << " not connections added" << std::endl;
        if (well_idx < -1) {
            std::cout << "Well not found " << wellname << std::endl;
        } else {
            std::cout << "Well " << wellname << " active " << wells_[well_idx].isActive()
                      << " fractures " << well_fractures_[well_idx].size() << std::endl;

            if (well_fractures_[well_idx].size() > 0) {
                for (const auto& frac : well_fractures_[well_idx]) {
                    std::cout << "Fracture " << frac.name() << " active " << frac.isActive()
                              << std::endl;
                }
            }
        }
    }

    return {};
}

template <typename Scalar>
void
FractureModel::assignGeomechWellState(WellState<Scalar>& wellState) const
{
    const auto nWells = this->wells_.size();
    for (auto i = 0 * nWells; i < nWells; ++i) {
        if (!wells_[i].isActive()) {
            continue;
        }

        const auto wsIx = wellState.index(this->wells_[i].name());
        if (!wsIx.has_value()) {
            continue;
        }

        auto& perfData = wellState[*wsIx].perf_data;

        if (perfData.connFracStatistics.size() != perfData.cell_index.size()) {
            perfData.connFracStatistics.resize(perfData.cell_index.size());
        }

        for (const auto& fracture : this->well_fractures_[i]) {
            if (!fracture.isActive()) {
                continue;
            }

            auto perfPos = std::find(
                perfData.cell_index.begin(), perfData.cell_index.end(), fracture.wellInfo().well_cell);
            if (perfPos == perfData.cell_index.end()) {
                continue;
            }

            // Possibly just "fracture.wellInfo().perf" instead.
            const auto perfIx = std::distance(perfData.cell_index.begin(), perfPos);
            fracture.assignGeomechWellState(perfData.connFracStatistics[perfIx]);
        }
    }
}

} // namespace Opm

// ===========================================================================
// Private member functions
// ===========================================================================

void
Opm::FractureModel::addFracturesPerpWell()
{
    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto& geom = elem.geometry();
            const auto& origin = geom.corner(1);
            const auto normal = origin - geom.corner(0);
            const auto elemIdx = emap.index(elem);

            this->well_fractures_[wellIx].emplace_back().init(
                fracWell.name(),
                /* connection index */ elemIdx,
                /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                /* global_index */ -1, // dummy for now
                /* well segment */ fracWell.segment(elemIdx),
                /* distance range */ fracWell.perfRange(elemIdx),
                /* seed origin */ origin,
                /* fracturing plane's normal vector */ normal,
                this->prm_);
        }
    }
}

void
Opm::FractureModel::addFracturesTensile()
{
    // https://link.springer.com/article/10.1007/s40948-023-00694-1

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            auto eigenVectors = Dune::FieldMatrix<double, 3, 3> {};
            auto eigenValues = Dune::FieldVector<double, 3> {};

            const auto elemIdx = emap.index(elem);

            const auto stressmat = ddm::symTensor2Matrix(fracWell.reservoirStress(elemIdx));
            Dune::FMatrixHelp::eigenValuesVectors(stressmat, eigenValues, eigenVectors);

            // Note: documentation for eigenValuesVectors() seems to imply
            // that minPos == eigenValues.begin() here.
            const auto minPos = std::min_element(eigenValues.begin(), eigenValues.end());
            const auto& normal = eigenVectors[std::distance(eigenValues.begin(), minPos)];

            this->well_fractures_[wellIx].emplace_back().init(
                fracWell.name(),
                /* connection index */ elemIdx,
                /* reservoir cell */ fracWell.reservoirCell(elemIdx),
                /* global_index */ -1, // dummy for now
                /* well segment */ fracWell.segment(elemIdx),
                /* distance range */ fracWell.perfRange(elemIdx),
                /* seed origin */ elem.geometry().corner(1),
                /* fracturing plane's normal vector */ normal,
                this->prm_);
        }
    }
}

namespace
{
std::unordered_map<int, std::size_t>
localSeedCells(const Opm::FractureWell& fracWell,
               const Opm::WellConnections& conns,
               const Opm::WellFractureSeeds& seeds)
{
    auto localSeedIxMap = std::unordered_map<int, std::size_t> {};

    auto connIx = [&fracWell, &conns](const std::size_t seedCellGlobal) {
        auto connPos = std::find_if(conns.begin(), conns.end(), [seedCellGlobal](const auto& conn) {
            return conn.global_index() == seedCellGlobal;
        });

        if (connPos == conns.end()) {
            return -1;
        }

        return fracWell.reservoirCell(std::distance(conns.begin(), connPos));
    };

    const auto& cells = seeds.seedCells();

    for (auto seedIx = 0 * cells.size(); seedIx < cells.size(); ++seedIx) {
        if (const auto ix = connIx(cells[seedIx]); ix >= 0) {
            localSeedIxMap.insert_or_assign(ix, seedIx);
        }
    }
    return localSeedIxMap;
}

} // Anonymous namespace

void
Opm::FractureModel::addFracturesWellSeed(const ScheduleState& sched)
{
    if (sched.wseed().empty()) {
        return;
    }

    using ElementMapper = Dune::MultipleCodimMultipleGeomTypeMapper<FractureWell::Grid::LeafGridView>;

    for (auto wellIx = 0 * this->wells_.size(); wellIx < this->wells_.size(); ++wellIx) {
        const auto& fracWell = this->wells_[wellIx];

        if (!sched.wseed.has(fracWell.name())) {
            continue;
        }

        const auto& wseed = sched.wseed(fracWell.name());
        const auto localSeeds
            = localSeedCells(fracWell, sched.wells(fracWell.name()).getConnections(), wseed);

        const auto emap = ElementMapper {fracWell.grid().leafGridView(), Dune::mcmgElementLayout()};

        for (const auto& elem : elements(fracWell.grid().leafGridView())) {
            const auto elemIdx = emap.index(elem);
            const auto seedPos = localSeeds.find(fracWell.reservoirCell(elemIdx));
            if (seedPos == localSeeds.end()) {
                continue;
            }

            // hack
            {
                const auto& seedSize = wseed.getSize(WellFractureSeeds::SeedIndex {seedPos->second});

                this->prm_.put("config.axis_scale", seedSize.verticalExtent());
                this->prm_.put("config.min_width", seedSize.width());
            }

            const auto& seedNormal = wseed.getNormal(WellFractureSeeds::SeedIndex {seedPos->second});

            const auto normal = Dune::FieldVector<double, 3> {
                seedNormal[0],
                seedNormal[1],
                seedNormal[2],
            };

            assert(normal.two_norm() > 0.0);

            const auto globalIndex = sched.wells(fracWell.name())
                .getConnections()[elemIdx].global_index();

            this->well_fractures_[wellIx].emplace_back()
                .init(fracWell.name(),
                      /* connection index */ elemIdx,
                      /* reservoir cell */ seedPos->first,
                      /* global_index */ static_cast<int>(globalIndex),
                      /* well segment */ fracWell.segment(elemIdx),
                      /* distance range */ fracWell.perfRange(elemIdx),
                      /* seed origin */ elem.geometry().corner(1),
                      /* fracturing plane's normal vector */ normal,
                      this->prm_);
        }
    }
}

// ===========================================================================
// Explicit specialisations.  No other code below separator.
// ===========================================================================

template void Opm::FractureModel::assignGeomechWellState(WellState<float>&) const;
template void Opm::FractureModel::assignGeomechWellState(WellState<double>&) const;
