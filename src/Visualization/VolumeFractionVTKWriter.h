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
 along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>
 */
/**
 * @file VolumeFractionVTKWriter.h
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief Debug output for volume fraction visualization
 **/

#ifndef VolumeFractionVTKWriter_h
#define VolumeFractionVTKWriter_h

#include <string>
#include <vector>
#include <list>
#include <map>

#include "../SubcellVFData.h"

namespace Clobscode
{
    class Quadrant;
    class MeshPoint;
    class Polyline;
    class QuadEdge;
    class EdgeInfo;

    class VolumeFractionVTKWriter
    {
    public:

        static bool writeQuadTreeWithVF(const std::string& name,
                                        const std::vector<Quadrant>& quadrants,
                                        const std::vector<MeshPoint>& points,
                                        const Polyline& geo);

        // Overload accepting any iterable container of Quadrant
        // (e.g. std::vector<Quadrant> or std::list<Quadrant> from the
        // TUSQH candidates list). Writes one file per call with the
        // s x s sample points of every quad as VTK_VERTEX cells, each
        // carrying its `winding_number` scalar (POINT_DATA) and the
        // parent q_id. This overload is defensive against cells whose
        // `mSampleSize` is set but whose `mWindingNumbers` was not
        // populated (e.g. cells produced by the legacy `-E`
        // IntersectionsVisitor path, which sets the state to Mixed
        // without recomputing the samples — those cells are skipped
        // silently).
        template <typename QuadrantContainer>
        static bool writeVFHeatmap(const std::string& name,
                                   const QuadrantContainer& quadrants,
                                   const std::vector<MeshPoint>& points);

        // Debug visualisation of the TUSQH winding-state classification:
        //   0 = Unknown
        //   1 = AllInside
        //   2 = AllOutside
        //   3 = Mixed
        // Plus provenance (QuadrantOrigin):
        //   0 = Classical
        //   1 = TusqhSplit
        //   2 = TusqhBalance
        //   3 = TusqhResolve
        //   4 = PreservedOutsider
        static bool writeWindingState(const std::string& name,
                                      const std::vector<Quadrant>& quadrants,
                                      const std::vector<MeshPoint>& points);

        // TUSQH §3.3 sub-cell volume fractions. Writes a single VTK
        // file containing, as POINT_DATA on the existing quadtree
        // vertices: subcell_vf, subcell_is_interior, subcell_sample_size.
        // As CELL_DATA on the quadtree quads: same metrics (the per-quad
        // subcell VF is the mean of the four corner vertices).
        static bool writeSubcellVertexVF(const std::string& name,
                                         const std::vector<Quadrant>& quadrants,
                                         const std::vector<MeshPoint>& points,
                                         double joinThreshold);

        // Writes the per-edge sub-cell volume fraction as a separate
        // UNSTRUCTURED_GRID of line cells (VTK_LINE = 3). Each line's
        // scalar is the edge's subcell VF.
        static bool writeSubcellEdgeVF(const std::string& name,
                                       const std::vector<Quadrant>& quadrants,
                                       const std::vector<MeshPoint>& points,
                                       const std::map<QuadEdge, EdgeSubcellVFData>& edgeSubcellVF,
                                       double joinThreshold);

        // Generic helper: writes the quadtree as a UNSTRUCTURED_GRID
        // and attaches a single CELL_DATA scalar field with the given
        // name and one double value per quad. Used for debug snapshots
        // that need to colour cells by an arbitrary integer label
        // (e.g. component id, keep/drop decision).
        //
        // `values.size()` must equal `quadrants.size()`.
        static bool writeQuadTreeWithCellArray(
            const std::string& name,
            const std::vector<Quadrant>& quadrants,
            const std::vector<MeshPoint>& points,
            const std::string& arrayName,
            const std::vector<double>& values);

        // Snapshot of the candidates list at a given TUSQH depth.
        // Writes the quadtree UNSTRUCTURED_GRID with:
        //   - winding_state  (per-cell, 0/1/2/3 as in writeWindingState)
        //   - volume_fraction (per-cell, NaN if not yet computed)
        //   - refinement_level (per-cell)
        //   - sample_size    (per-cell)
        // The name suffix is appended verbatim (typically "_tusqh_iter<N>"
        // for the main loop, "_tusqh_resolve_iter<N>" for the resolve
        // pass).
        template <typename QuadrantContainer>
        static bool writeCandidatesSnapshot(const std::string& name,
                                             const QuadrantContainer& candidates,
                                             const std::vector<MeshPoint>& points,
                                             const std::string& suffix);
    };
}

#endif