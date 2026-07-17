/*
  <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

  Copyright (C) <2013,2026>  <Claudio Lobos, Fabrice Jaillet, Felipe Marchant>
  All rights reserved.

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
* @file PinchVisitor.h
* @brief Abstract base class for visiting quadtree vertices (paper §3.4).
*
* Unlike the existing `Visitor` hierarchy which operates on `Quadrant`,
* pinch detection operates on the vertices of the cubical complex.
* Each vertex is visited once with the list of incident quads, so the
* concrete visitor can classify the local 2x2 (or hanging) neighbourhood.
**/

#ifndef PinchVisitor_h
#define PinchVisitor_h 1

#include <vector>
#include <map>

namespace Clobscode
{
    // Forward declarations to avoid heavy includes.
    class MeshPoint;
    class Quadrant;
    class QuadEdge;
    struct EdgeInfo;

    //--------------------------------------------------------------------
    // PinchVisitor
    //
    // Standalone visitor (does NOT extend the existing `Visitor` base
    // because pinch operates on vertices, not quadrants). The concrete
    // subclass implements `visit(vertexId, incidentQuads)` to inspect
    // the local neighbourhood and update the per-vertex pinch state.
    //
    // Setters are provided for the read-only context the visitor needs:
    //   - points:     the mesh points (for vertex coordinates and
    //                 sub-cell VF at the vertex).
    //   - quadrants:  the cubical complex (to read incidentQuads).
    //   - mapEdges:   the edge map (for the consistency graph in Fase 4).
    //   - threshold:  the join threshold used by the resolution step
    //                 to decide Connect vs Separate.
    //
    // `visit` returns true if iteration should continue. Returning false
    // stops the iteration early (useful for unit tests).
    //--------------------------------------------------------------------
    class PinchVisitor
    {
    public:
        PinchVisitor();
        virtual ~PinchVisitor();

        void setPoints(const std::vector<MeshPoint> *mp);
        void setQuadrants(const std::vector<Quadrant> *q);
        void setMapEdges(const std::map<QuadEdge, EdgeInfo> *me);
        void setJoinThreshold(double t);

        // Called for each vertex of the quadtree. `vertexId` is the
        // index into the `points` vector. `incidentQuads` is the list
        // of Quadrant indices whose pointIndex list contains `vertexId`.
        //
        // Returns true to continue iteration, false to stop early.
        virtual bool visit(unsigned int vertexId,
                           const std::vector<unsigned int> &incidentQuads) = 0;

    protected:
        const std::vector<MeshPoint> *mPoints;
        const std::vector<Quadrant> *mQuadrants;
        const std::map<QuadEdge, EdgeInfo> *mMapEdges;
        double mJoinThreshold;
    };

} // namespace Clobscode

#endif // PinchVisitor_h
