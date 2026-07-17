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
* @file PinchCase.h
* @brief Classification of grid vertices regarding pinch configurations
*        (paper TUSQH §3.4, Fig. 9).
*
* This file defines the data structures used by the pinch detection
* and resolution pipeline. The pinch detection step (Fase 2) writes
* a `PinchCase` per vertex; the resolution step (Fase 3-5) consumes
* it to decide whether to apply the growing or the shrinking template.
**/

#ifndef PinchCase_h
#define PinchCase_h 1

#include <vector>

namespace Clobscode
{

    //--------------------------------------------------------------------
    // PinchCase
    //
    // Classification of a single grid vertex regarding its neighbourhood
    // in the cubical complex. Mirrors paper Fig. 9.
    //
    //   None               - not a pinch (default).
    //   Vertex_2x2_Chess   - 2D vertex pinch (Fig. 9a): exactly 4 incident
    //                        quads, 2 diagonal interior + 2 diagonal exterior.
    //                        The two interior quads meet at this vertex
    //                        without sharing an edge (chess-pattern).
    //   Vertex_Hanging_2   - 2 incident quads only (T-junction, Fig. 9c/d).
    //                        Either same level (edge) or different levels
    //                        (asymmetric pinch, paper §3.4 second paragraph).
    //   Vertex_Hanging_3   - 3 incident quads (T-junction corner, Fig. 9e).
    //   Vertex_Corner_1    - 1 incident quad (domain corner, Fig. 9l).
    //                        Cannot be a pinch (always a manifold boundary).
    //--------------------------------------------------------------------
    enum class PinchCase {
        None              = 0,
        Vertex_2x2_Chess  = 1,
        Vertex_Hanging_2  = 2,
        Vertex_Hanging_3  = 3,
        Vertex_Corner_1   = 4
    };

    //--------------------------------------------------------------------
    // PinchResolution
    //
    // Which template to apply when repairing a pinch.
    //   None       - no action (no pinch, or already resolved).
    //   Connect    - paper Fig. 12 (growing template): the pinch vertex is
    //                INTERIOR according to sub-cell VF; bridge the two
    //                diagonal components with a 1-to-5 split and KEEP the
    //                bridge quad.
    //   Separate   - paper Fig. 11 (shrinking template): the pinch vertex
    //                is EXTERIOR according to sub-cell VF; bridge the two
    //                diagonal components with a 1-to-5 split and REMOVE
    //                the child quad that contains the pinch vertex, so the
    //                two components end up disconnected.
    //--------------------------------------------------------------------
    enum class PinchResolution {
        None      = 0,
        Connect   = 1,
        Separate  = 2
    };

    //--------------------------------------------------------------------
    // VertexPinchInfo
    //
    // Per-vertex state produced by `PinchDetector::detectAll` and consumed
    // by the resolution passes (Fases 3, 4, 5).
    //
    // Fields:
    //   pinCase           - which kind of pinch (or None).
    //   resolution        - which template to apply (Connect / Separate).
    //                       Filled by Fase 2 from sub-cell VF at the vertex.
    //   incidentQuads     - indices into the Quadrants vector of the
    //                       quads that contain this vertex. 1 to 4 entries.
    //   interiorDiagA,B   - for Vertex_2x2_Chess: indices of the two
    //                       interior diagonal quads. The two exterior
    //                       diagonal quads are at positions
    //                       incidentQuads \ {interiorDiagA, interiorDiagB}.
    //   subcellVF         - sub-cell volume fraction at this vertex
    //                       (from mVertexSubcellVF), used to choose
    //                       Connect vs Separate in auto-resolution mode.
    //   equalizedLevel    - for Vertex_Hanging_2 at different levels:
    //                       refinement level that the larger quad must be
    //                       subdivided to in order to apply the standard
    //                       1-to-5 template (paper §3.4, Fig. 9d asymmetric
    //                       case, "pre-select" / equalize step).
    //--------------------------------------------------------------------
    struct VertexPinchInfo {
        PinchCase pinCase = PinchCase::None;
        PinchResolution resolution = PinchResolution::None;
        std::vector<unsigned int> incidentQuads;
        unsigned int interiorDiagA = 0;
        unsigned int interiorDiagB = 0;
        double subcellVF = 0.0;
        unsigned short equalizedLevel = 0;
    };

} // namespace Clobscode

#endif // PinchCase_h
