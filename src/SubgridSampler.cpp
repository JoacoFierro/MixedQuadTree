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
  * @file SubgridSampler.cpp
  * @brief See SubgridSampler.h.
  **/

#include "SubgridSampler.h"
#include "Point3D.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    unsigned int SubgridSampler::sanitizeSampleSize(unsigned int s)
    {
        if (s < 2) return 2;
        if (s % 2 != 0) return s + 1;
        return s;
    }

    //--------------------------------------------------------------------------------
    SubgridSampler::Result SubgridSampler::sampleVertex(
        const Point3D& vertexPos,
        const vector<pair<unsigned int, double>>& incidentEdges,
        const Polyline& input,
        unsigned int s) const
    {
        Result r;
        r.sampleSize = sanitizeSampleSize(s);

        // Fictitious cell side = max length of incident edges. If the
        // vertex has no incident edges (shouldn't happen in a quadtree
        // interior point) we fall back to a tiny epsilon to avoid a
        // zero-size cell that would produce NaN samples.
        double side = 0.0;
        for (const auto& e : incidentEdges) {
            if (e.second > side) side = e.second;
        }
        if (side <= 0.0) side = std::numeric_limits<double>::epsilon();
        r.cellSizeAlong = side;
        r.cellSizePerpendicular = side;

        const double halfSide = 0.5 * side;
        const double step = side / static_cast<double>(r.sampleSize);

        r.windingNumbers.reserve(r.sampleSize * r.sampleSize);

        for (unsigned int i = 0; i < r.sampleSize; ++i) {
            for (unsigned int j = 0; j < r.sampleSize; ++j) {
                // Sample at the centre of each sub-cell of the fictitious
                // grid, exactly like Quadrant::getSamplePoint does for
                // maximal-dimension cells (so the same convention is used
                // for all dimensions).
                double x = vertexPos[0] - halfSide + (i + 0.5) * step;
                double y = vertexPos[1] - halfSide + (j + 0.5) * step;
                Point3D sample(x, y, 0.0);
                int wn = input.windingNumber(sample);
                r.windingNumbers.push_back(static_cast<double>(wn));
            }
        }

        double sum = 0.0;
        for (double v : r.windingNumbers) sum += v;
        r.volumeFraction = r.windingNumbers.empty()
                              ? 0.0
                              : sum / static_cast<double>(r.windingNumbers.size());
        return r;
    }

    //--------------------------------------------------------------------------------
    SubgridSampler::Result SubgridSampler::sampleEdge(
        const Point3D& edgeA,
        const Point3D& edgeB,
        const vector<double>& quadPerpThickness,
        const Polyline& input,
        unsigned int s) const
    {
        Result r;
        r.sampleSize = sanitizeSampleSize(s);

        // Edge length and midpoint.
        double dx = edgeB[0] - edgeA[0];
        double dy = edgeB[1] - edgeA[1];
        double edgeLen = std::sqrt(dx * dx + dy * dy);
        if (edgeLen <= 0.0) edgeLen = std::numeric_limits<double>::epsilon();
        r.cellSizeAlong = edgeLen;

        // Perpendicular thickness = max over both sides. If neither side
        // reports a quad (boundary edge) we still keep the cell finite by
        // using the edge length as a sensible default.
        double perp = 0.0;
        for (double t : quadPerpThickness) {
            if (t > perp) perp = t;
        }
        if (perp <= 0.0) perp = edgeLen;
        r.cellSizePerpendicular = perp;

        // Build an orthonormal frame: along = edge direction, perp = rotate 90°.
        double ax = dx / edgeLen;
        double ay = dy / edgeLen;
        double px = -ay;
        double py =  ax;

        double mx = 0.5 * (edgeA[0] + edgeB[0]);
        double my = 0.5 * (edgeA[1] + edgeB[1]);

        const double halfAlong = 0.5 * edgeLen;
        const double halfPerp  = 0.5 * perp;
        const double stepAlong = edgeLen   / static_cast<double>(r.sampleSize);
        const double stepPerp  = perp      / static_cast<double>(r.sampleSize);

        r.windingNumbers.reserve(r.sampleSize * r.sampleSize);

        for (unsigned int i = 0; i < r.sampleSize; ++i) {
            for (unsigned int j = 0; j < r.sampleSize; ++j) {
                double u = -halfAlong + (i + 0.5) * stepAlong;
                double v = -halfPerp  + (j + 0.5) * stepPerp;
                double x = mx + u * ax + v * px;
                double y = my + u * ay + v * py;
                Point3D sample(x, y, 0.0);
                int wn = input.windingNumber(sample);
                r.windingNumbers.push_back(static_cast<double>(wn));
            }
        }

        double sum = 0.0;
        for (double v : r.windingNumbers) sum += v;
        r.volumeFraction = r.windingNumbers.empty()
                              ? 0.0
                              : sum / static_cast<double>(r.windingNumbers.size());
        return r;
    }

    //--------------------------------------------------------------------------------
    vector<pair<unsigned int, double>>
    SubgridSampler::buildIncidentEdgeList(unsigned int vertexIdx,
                                          const map<QuadEdge, EdgeInfo>& mapEdges,
                                          const vector<MeshPoint>& points)
    {
        vector<pair<unsigned int, double>> result;
        const Point3D& p = points[vertexIdx].getPoint();

        for (const auto& entry : mapEdges) {
            const QuadEdge& e = entry.first;
            const unsigned int a = e[0];
            const unsigned int b = e[1];
            if (a != vertexIdx && b != vertexIdx) continue;

            unsigned int other = (a == vertexIdx) ? b : a;
            const Point3D& q = points[other].getPoint();
            double len = (p - q).Norm();
            result.emplace_back(other, len);
        }
        return result;
    }

    //--------------------------------------------------------------------------------
    vector<double>
    SubgridSampler::buildQuadPerpThickness(const QuadEdge& edge,
                                           const map<QuadEdge, EdgeInfo>& mapEdges,
                                           const unordered_map<unsigned int, unsigned int>& qIdToIdx,
                                           const vector<Quadrant>& quadrants,
                                           const vector<MeshPoint>& points)
    {
        vector<double> result;
        auto it = mapEdges.find(edge);
        if (it == mapEdges.end()) return result;

        const EdgeInfo& info = it->second;
        // EdgeInfo stores (midpointIdx, q1, q2) where q1/q2 are q_ids
        // (not vector indices). The q_id is set at construction time
        // by SplitVisitor / Quadrant ctor and never reassigned, while
        // the Quadrants vector may be reordered during the pipeline
        // (e.g. after `resolveArchipelagos`'s compact-and-rebuild, or
        // even during windingSubdivide's iteration over `idx_pos_map`).
        // We must look up the vector position via qIdToIdx before
        // indexing `quadrants[]`. See BUGS_FOUND.md Issue #9.
        for (unsigned int k = 1; k <= 2; ++k) {
            unsigned int qId = info[k];
            if (qId == std::numeric_limits<unsigned int>::max()) continue;
            auto qIt = qIdToIdx.find(qId);
            if (qIt == qIdToIdx.end()) continue;
            const Quadrant& q = quadrants[qIt->second];

            // The quadrant's parent corners are at indices 0..3 of
            // getPointIndex(). Sub-elements may further split it. We
            // approximate the perpendicular thickness by the distance
            // from the edge midpoint to the centroid of the quadrant.
            // For a convex quad this equals half the perpendicular
            // extent up to O(h) error, which is good enough for VF.
            const Point3D& A = points[edge[0]].getPoint();
            const Point3D& B = points[edge[1]].getPoint();
            double mx = 0.5 * (A[0] + B[0]);
            double my = 0.5 * (A[1] + B[1]);
            Point3D edgeMid(mx, my, 0.0);

            Point3D c(0,0,0);
            const auto& pi = q.getPointIndex();
            if (pi.empty()) continue;
            for (unsigned int v : pi) c += points[v].getPoint();
            c /= static_cast<double>(pi.size());

            double proj = (c - edgeMid).Norm();
            if (proj > 0.0) result.push_back(proj);
        }
        return result;
    }
}
