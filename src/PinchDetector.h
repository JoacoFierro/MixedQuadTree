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
* @file PinchDetector.h
* @brief Detects pinch cases (paper TUSQH §3.4 / Fig. 9) on the
*        cubical complex after TUSQH subdivision and bridge joining.
*
* Fase 1 introduces the API and stub implementation. Fase 2 fills in
* the body of `detectAtVertex` for the 2x2 chess pattern (Fig. 9a).
* Fase 5 adds the hanging-node branches.
**/

#ifndef PinchDetector_h
#define PinchDetector_h 1

#include <vector>
#include <map>
#include "PinchCase.h"

namespace Clobscode
{
    // Forward declarations.
    class MeshPoint;
    class Quadrant;
    class Polyline;
    class QuadEdge;
    struct EdgeInfo;

    //--------------------------------------------------------------------
    // PinchDetector
    //
    // Detects pinch configurations on the quadtree and stores the result
    // in a `VertexPinchInfo` vector indexed by vertex id.
    //
    // Usage:
    //   PinchDetector pd;
    //   pd.setPoints(&points);
    //   pd.setQuadrants(&quadrants);
    //   pd.setMapEdges(&MapEdges);
    //   pd.setJoinThreshold(joinThreshold);
    //   unsigned int n = pd.detectAll();
    //   const auto &info = pd.getPinchInfo();
    //
    // After `detectAll`:
    //   - getPinchInfo()[i]  - state of vertex i.
    //   - getPinchCount()    - number of vertices whose pinCase != None.
    //
    // The detector does NOT mutate the cubical complex. The resolution
    // passes (Fases 3, 4, 5) read the produced info and apply templates.
    //--------------------------------------------------------------------
    class PinchDetector
    {
    public:
        PinchDetector();

        void setPoints(const std::vector<MeshPoint> *mp);
        void setQuadrants(const std::vector<Quadrant> *q);
        void setMapEdges(const std::map<QuadEdge, EdgeInfo> *me);
        void setPolyline(const Polyline *input);
        void setJoinThreshold(double t);
        void setSampleSize(unsigned int s);

        // Iterates over every vertex of the quadtree and classifies it.
        // Returns the number of vertices whose pinCase != None.
        unsigned int detectAll();

        const std::vector<VertexPinchInfo> &getPinchInfo() const
        {
            return mPinchInfo;
        }

        unsigned int getPinchCount() const;

    private:
        // Finds the indices of the Quadrants whose pointIndex list
        // contains `vertexId`. Returns 1, 2, 3 or 4 entries.
        std::vector<unsigned int> findIncidentQuads(unsigned int vertexId) const;

        // Classifies a single vertex. Implementation in Fase 2.
        PinchCase detectAtVertex(unsigned int vertexId,
                                 const std::vector<unsigned int> &incidentQuads);

        // Assigns `resolution` (Connect/Separate) based on sub-cell VF
        // at the vertex. Called by `detectAll` after `detectAtVertex`.
        void assignResolution(unsigned int vertexId);

        const std::vector<MeshPoint> *mPoints;
        const std::vector<Quadrant> *mQuadrants;
        const std::map<QuadEdge, EdgeInfo> *mMapEdges;
        const Polyline *mPolyline;
        double mJoinThreshold;
        unsigned int mSampleSize;

        std::vector<VertexPinchInfo> mPinchInfo;
    };

} // namespace Clobscode

#endif // PinchDetector_h
