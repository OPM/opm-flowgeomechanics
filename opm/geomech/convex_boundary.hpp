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

#ifndef OPM_GEOMECH_CONVEX_BOUNDARY_HPP_INCLUDED
#define OPM_GEOMECH_CONVEX_BOUNDARY_HPP_INCLUDED

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <vector>

namespace Opm
{
// ----------------------------------------------------------------------------
// @brief Template function to find the best 2D plane for projecting a set of 3D
// points onto.
//
// @tparam Point3D A type representing a 3-dimensional point. It is expected to
// support indexing via the subscript operator (e.g., pts[i][d]).
//
// @param pts A constant reference to a vector of 3D points. The vector must
// contain at least one point.
//
// @return std::array<int, 2> An array containing the indices of the two
// dimensions that form the best 2D plane.
template <typename Point3D>
std::array<int, 2>
best2Dplane(const std::vector<Point3D>& pts)
// ----------------------------------------------------------------------------
{
    assert(!pts.empty());

    const std::size_t N = pts.size();
    std::array<double, 3> min_coords {pts[0][0], pts[0][1], pts[0][2]};
    std::array<double, 3> max_coords {pts[0][0], pts[0][1], pts[0][2]};

    for (std::size_t i = 0; i != N; ++i) {
        for (std::size_t d = 0; d != 3; ++d) {
            min_coords[d] = std::min(min_coords[d], pts[i][d]);
            max_coords[d] = std::max(max_coords[d], pts[i][d]);
        }
    }

    for (std::size_t d = 0; d != 3; ++d) {
        max_coords[d] -= min_coords[d];
    }

    const std::size_t eliminated_dim
        = std::min_element(max_coords.begin(), max_coords.end()) - max_coords.begin();

    return eliminated_dim == 0 ? std::array<int, 2> {1, 2}
        : eliminated_dim == 1  ? std::array<int, 2> {0, 2}
                               : std::array<int, 2> {0, 1};
}

// ----------------------------------------------------------------------------
// This function projects a set of 3D points onto a specified 2D plane. The
// plane is defined by two dimensions (indices) from the 3D space.
template <typename Point3D>
std::vector<double>
projectPointsTo2DPlane(const std::vector<Point3D>& pts, const std::array<int, 2>& plane)
// ----------------------------------------------------------------------------
{
    std::vector<double> result;

    for (const auto& p : pts) {
        result.push_back(p[plane[0]]);
        result.push_back(p[plane[1]]);
    }

    return result;
}

// ----------------------------------------------------------------------------
void
search_cv_points(const std::size_t& start,
                 const std::size_t& end,
                 const std::vector<double>& pts2D,
                 std::vector<int>& flag,
                 const std::vector<std::size_t>& candidates)
// ----------------------------------------------------------------------------
{
    const double x0 = pts2D[2 * start], y0 = pts2D[2 * start + 1];
    const double x1 = pts2D[2 * end], y1 = pts2D[2 * end + 1];
    const double dx0 = x1 - x0, dy0 = y1 - y0;
    const double norm = std::sqrt(dx0 * dx0 + dy0 * dy0);
    const double nx = dx0 / norm, ny = dy0 / norm;

    double max_dist = 0;

    // identify all points to the right of the line defined by start and end,
    // as well as the candidate furthest from the line.
    std::vector<std::size_t> right_side_candidates;
    std::size_t best_ix = 0;
    for (const auto& ix : candidates) {
        if (ix == start || ix == end) {
            continue;
        }

        const double x = pts2D[2 * ix], y = pts2D[2 * ix + 1];
        const double dx = x - x0, dy = y - y0;
        const double dist = dx * ny - dy * nx;

        if (dist > 0) {
            right_side_candidates.push_back(ix);

            if (dist > max_dist) {
                max_dist = dist;
                best_ix = ix;
            }
        }
    }

    if (!right_side_candidates.empty()) {
        // insert into convex hull sequence
        flag[start] = best_ix;
        flag[best_ix] = end;

        search_cv_points(start, best_ix, pts2D, flag, right_side_candidates);
        search_cv_points(best_ix, end, pts2D, flag, right_side_candidates);
    }
}

// ----------------------------------------------------------------------------
// @brief This function finds the convex hull of a set of 2D points.
// It returns the indices of the points that form the convex hull.
std::vector<std::size_t>
convex_hull(const std::vector<double>& pts2D)
// ----------------------------------------------------------------------------
{
    assert(!pts2D.empty());

    const std::size_t N = pts2D.size() / 2; // number of 2D points

    // identify points touching the bounding box - these must be on the convex hull
    std::array<std::size_t, 4> init_pts {0, 0, 0, 0}; // minx, miny, max, maxy
    std::array<double, 2> xbounds {pts2D[0], pts2D[0]}, ybounds {pts2D[1], pts2D[1]};

    for (std::size_t ix = 1; ix != N; ++ix) {
        if (pts2D[2 * ix] < xbounds[0]) {
            xbounds[0] = pts2D[2 * ix];
            init_pts[0] = ix;
        }

        if (pts2D[2 * ix] > xbounds[1]) {
            xbounds[1] = pts2D[2 * ix];
            init_pts[2] = ix;
        }

        if (pts2D[2 * ix + 1] < ybounds[0]) {
            ybounds[0] = pts2D[2 * ix + 1];
            init_pts[1] = ix;
        }

        if (pts2D[2 * ix + 1] > ybounds[1]) {
            ybounds[1] = pts2D[2 * ix + 1];
            init_pts[3] = ix;
        }
    }

    std::vector<int> flag(N, -1);
    for (int i = 0; i != 4; ++i) {
        flag[init_pts[i]] = init_pts[(i + 1) % 4]; // each cv point points to the next found
    }

    std::vector<std::size_t> candidates(N, 0);
    std::iota(candidates.begin(), candidates.end(), 0);

    const std::size_t start_cv = init_pts[0];
    bool at_start = true;
    for (std::size_t run_cv = start_cv; at_start || run_cv != start_cv;
         run_cv = flag[run_cv], at_start = false) {
        search_cv_points(run_cv, flag[run_cv], pts2D, flag, candidates);
    }

    // collect all indentified convex hull points, in order
    std::vector<std::size_t> result(1, start_cv);
    for (std::size_t next_cv = flag[start_cv]; next_cv != start_cv; next_cv = flag[next_cv]) {
        result.push_back(next_cv);
    }

    return result;
}

} // end namespace Opm

#endif // OPM_GEOMECH_CONVEX_BOUNDARY_HPP_INCLUDED
