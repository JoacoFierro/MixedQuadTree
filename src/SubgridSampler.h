/*
  <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

  Copyright (C) <2013,2024>  <Claudio Lobos, Fabrice Jaillet> All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Lesser General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>.
*/
/**
  * @file SubgridSampler.h
  * @author Felipe Marchant
  * @brief Implements the "subgrid sampling" technique from the TUSQH paper
  *        (Shawcroft, Shepherd & Mitchell, IMR 2025) for cells of dimension
  *        lower than the maximal one (vertices and edges in 2D).
  *
  *        For an n-cell we define a fictitious cell centered on it. The
  *        volume fraction of that n-cell is the mean winding number of an
  *        s x s grid of sample points uniformly distributed inside the
  *        fictitious cell. We use only EVEN sample sizes (s>=2) so the
  *        samples never lie on the cell boundary (paper Fig. 8).
  *
  *        0-cell (vertex): fictitious cell is the square centered at the
  *        vertex with side = max length of incident edges.
  *
  *        1-cell (edge): fictitious cell is the rectangle centered on the
  *        edge midpoint. Length along the edge = the edge length.
  *        Perpendicular extent = max perpendicular thickness of the
  *        adjacent quadrants (0 if the edge has only one adjacent quad,
  *        i.e. it sits on the domain boundary).
  *
  *        This class is stateless and side-effect free: each method just
  *        computes the requested volume fraction and returns it.
  **/

#ifndef SubgridSampler_h
#define SubgridSampler_h 1

#include "MeshPoint.h"
#include "Polyline.h"
#include "QuadEdge.h"
#include "EdgeInfo.h"
#include "Quadrant.h"

#include <map>
#include <unordered_map>
#include <vector>
#include <utility>

using std::map;
using std::vector;
using std::pair;

namespace Clobscode
{
    class SubgridSampler
    {
    public:

        /// Result of a subgrid sampling: the mean wn, the raw wn values
        /// (s*s of them), and a few diagnostics (sample size + cell size).
        struct Result {
            double volumeFraction;
            vector<double> windingNumbers;
            unsigned int sampleSize;
            // Side length of the fictitious cell used for sampling.
            // For an edge this is {lengthAlong, lengthPerpendicular};
            // for a vertex this is {side, side}.
            double cellSizeAlong;
            double cellSizePerpendicular;
        };

        SubgridSampler() = default;
        virtual ~SubgridSampler() = default;

        /// s must be even and >= 2 (paper restriction). Odd values are
        /// promoted to the next even integer.
        static unsigned int sanitizeSampleSize(unsigned int s);

        /// Compute the subcell volume fraction for a vertex.
        /// `incidentEdges` is the list of (other_endpoint_idx, length) for
        /// every quadtree edge that touches the vertex. The fictitious
        /// cell size is max(length) of those edges (0 if no edges).
        Result sampleVertex(const Point3D& vertexPos,
                            const vector<pair<unsigned int, double>>& incidentEdges,
                            const Polyline& input,
                            unsigned int s) const;

        /// Compute the subcell volume fraction for an edge.
        /// `quadPerpThickness` is the perpendicular thickness of each
        /// adjacent quadrant (0 if no quad, e.g. boundary edge).
        Result sampleEdge(const Point3D& edgeA,
                          const Point3D& edgeB,
                          const vector<double>& quadPerpThickness,
                          const Polyline& input,
                          unsigned int s) const;

        // ----- Convenience helpers used by Mesher -----

        /// Build the list of incident edge lengths for a given vertex by
        /// scanning the quadtree edge map. Returns empty list if the
        /// vertex has no incident edges.
        static vector<pair<unsigned int, double>>
        buildIncidentEdgeList(unsigned int vertexIdx,
                             const map<QuadEdge, EdgeInfo>& mapEdges,
                             const vector<MeshPoint>& points);

        /// For an edge, return the perpendicular thickness of every
        /// adjacent quadrant (0 if a quadrant does not exist or if the
        /// quadrilateral is degenerate). The returned vector has up to 2
        /// entries (one per side).
        ///
        /// `qIdToIdx` maps each quadrant's q_id (== EdgeInfo info[1/2])
        /// to its position in the `quadrants` vector. EdgeInfo stores
        /// q_ids (not vector indices), so the lookup is required.
        /// Build it once before the edge loop:
        ///   unordered_map<unsigned int, unsigned int> qIdToIdx;
        ///   for (unsigned int i = 0; i < quadrants.size(); ++i)
        ///       qIdToIdx[quadrants[i].getIndex()] = i;
        static vector<double>
        buildQuadPerpThickness(const QuadEdge& edge,
                               const map<QuadEdge, EdgeInfo>& mapEdges,
                               const unordered_map<unsigned int, unsigned int>& qIdToIdx,
                               const vector<Quadrant>& quadrants,
                               const vector<MeshPoint>& points);
    };
}

#endif
