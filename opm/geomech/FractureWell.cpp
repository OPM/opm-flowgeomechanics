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

#include <opm/geomech/FractureWell.hpp>

#include <dune/common/fvector.hh>

#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/foamgrid/foamgrid.hh>

#include <opm/models/io/vtkmultiwriter.hh>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace Opm
{
FractureWell::FractureWell(const std::string& outputprefix,
                           const std::string& casename,
                           const std::string& name,
                           const std::vector<Connection>& conns,
                           const std::vector<Point3D>& points,
                           const std::vector<Segment>& segments)
    : outputprefix_ {outputprefix}
    , casename_ {casename}
    , name_ {name}
    , conns_ {conns}
{
    this->init(points, segments);

    this->vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(),
                                                                             Dune::VTK::nonconforming);

    reservoir_stress_.resize(this->conns_.size());
    std::fill(reservoir_stress_.begin(), reservoir_stress_.end(),
              100e5); // random initialization
}

void
FractureWell::init(const std::vector<Point3D>& points, const std::vector<Segment>& segments)
{
    Dune::GridFactory<Grid> factory;

    for (std::size_t i = 0; i < points.size(); ++i) {
        factory.insertVertex(points[i]);
    }

    for (std::size_t i = 0; i < segments.size(); ++i) {
        std::vector<unsigned int> seg {segments[i][0], segments[i][1]};

        factory.insertElement(Dune::GeometryTypes::line, seg);
    }

    grid_ = factory.createGrid();

    this->resetWriters();
}

void
FractureWell::write() const
{
    const std::string filename = outputprefix_ + "/" + this->name();
    vtkwriter_->write(filename.c_str());
}

void
FractureWell::resetWriters()
{
    // nead to be reseat if grid is changed ??
    vtkwriter_ = std::make_unique<Dune::VTKWriter<Grid::LeafGridView>>(grid_->leafGridView(),
                                                                       Dune::VTK::nonconforming);

    const std::string outputdir = outputprefix_;
    const std::string simName = casename_ + "_" + this->name();
    const std::string multiFileName = "";

    vtkmultiwriter_ = std::make_unique<Opm::VtkMultiWriter<Grid::LeafGridView, VTKFormat>>(
        /*async*/ false, grid_->leafGridView(), outputdir, simName, multiFileName);
}

void
FractureWell::writemulti(double time) const
{
    std::vector<double> reservoir_pressure = reservoir_pressure_;
    std::vector<double> reservoir_temperature = reservoir_temperature_;
    std::vector<double> perf_pressure = perf_pressure_;

    vtkmultiwriter_->beginWrite(time);

    if (!perf_pressure.empty()) {
        vtkmultiwriter_->attachScalarElementData(perf_pressure, "PerfPressure");
    }

    if (!reservoir_pressure.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_pressure, "ReservoirPressure");
    }

    if (!reservoir_temperature.empty()) {
        vtkmultiwriter_->attachScalarElementData(reservoir_temperature, "ReservoirTemperature");
    }

    vtkmultiwriter_->endWrite(false);
}

} // namespace Opm
