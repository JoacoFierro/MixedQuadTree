/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2020>  <Claudio Lobos> All rights reserved.

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
* @file Mesher.h
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#ifndef Mesher_h
#define Mesher_h 1

#include "Polyline.h"
#include "FEMesh.h"
#include "GridMesher.h"
#include "Quadrant.h"
#include "QuadEdge.h"
#include "EdgeInfo.h"
#include "Services.h"
#include "RefinementRegion.h"
#include "RefinementCubeRegion.h"

#include "Visitors/SplitVisitor.h"
#include "Visitors/IntersectionsVisitor.h"
#include "Visitors/TransitionPatternVisitor.h"
#include "Visitors/SurfaceTemplatesVisitor.h"
#include "Visitors/RemoveSubElementsVisitor.h"
#include "Visitors/EdgeVisitor.h"
#include "Visitors/OneIrregularVisitor.h"
#include "Visitors/WindingNumberVisitor.h"
#include "Visitors/WindingNumberSubdivisionVisitor.h"
#include "Visualization/VolumeFractionVTKWriter.h"
#include "SubgridSampler.h"
#include "SubcellVFData.h"
// ----- Modificado por Joaquin Fierro --------------
#include "Visitors/QuadAliasing.h"
// ----- Fin modificacion ---------------------------

#include <list>
#include <vector>
#include <set>
#include <map>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <chrono>

using std::vector;
using std::list;
using std::set;
using Clobscode::QuadEdge;
using Clobscode::Polyline;
using Clobscode::RefinementRegion;

#define VTKOUT true //CL Debbuging

namespace Clobscode
{
	
	class Mesher{
		
	public:

		Mesher();
		
		virtual ~Mesher();
				
        virtual std::shared_ptr<FEMesh> generateMesh(Polyline &input, const unsigned short &rl,
                                  const string &name, list<RefinementRegion *> &all_reg,
                                  const bool &debugging, unsigned int sampleSize=2,
                                  bool decoration=false,
                                  bool useTusqh=false,
                                  unsigned int tusqhSampleSize=2,
                                  bool refineOnEdgeIntersect=false,
                                  unsigned int tusqhExtraResolveDepth=0,bool Aliasing=false,
                                  bool useSubgrid=false,
                                  unsigned int subgridSampleSize=2,
                                  double subgridJoinThreshold=0.5,
                                  unsigned int subgridMinComponentCells=5);

        virtual std::shared_ptr<FEMesh> refineMesh(Polyline &input, const unsigned short &rl,
                                  const string &name, list<unsigned int> &roctli,
                                  list<RefinementRegion *> &all_reg,
                                  GeometricTransform &gt, const unsigned short &minrl,
                                  const unsigned short &omaxrl, const bool &debugging,
                                  unsigned int sampleSize=2, bool decoration=false,
                                  bool useTusqh=false,
                                  unsigned int tusqhSampleSize=2,
                                  bool refineOnEdgeIntersect=false,
                                  unsigned int tusqhExtraResolveDepth=0,bool Aliasing=false,
                                  bool useSubgrid=false,
                                  unsigned int subgridSampleSize=2,
                                  double subgridJoinThreshold=0.5,
                                  unsigned int subgridMinComponentCells=5);

        
        virtual void setInitialState(vector<MeshPoint> &epts, vector<Quadrant> &eocts,
                                     map<QuadEdge, EdgeInfo> &medgs);
        
	protected:
        
        virtual void splitQuadrants(const unsigned short &rl, Polyline &input,
                                    list<unsigned int> &roctli,
                                    list<RefinementRegion *> &all_reg, const string &name,
                                    const unsigned short &minrl, const unsigned short &maxrl,
                                    const bool &debugging=false);

        virtual void generateQuadtreeMesh(const unsigned short &rl, Polyline &input,
                                          const list<RefinementRegion *> &all_reg,
                                          const string &name, const unsigned short &minrl,
                                          const unsigned short &givenmaxrl=0,
                                          const bool &debugging=false,
                                          unsigned int new_q_idx=0,bool Aliasing=false);

        // TUSQH subdivision loop. Subdivides each quadrant while its s x s
        // sample-point winding numbers are ambiguous (i.e. some samples are
        // inside and some outside the Polyline). Max depth is bounded by
        // `maxDepth`.
        // Optional legacy mode: when refineOnEdgeIntersect is true the
        // cell is also refined whenever it geometrically intersects the
        // input Polyline, regardless of the winding criterion.
        virtual void windingSubdivide(Polyline &input, unsigned int maxDepth,
                                      unsigned int tusqhSampleSize,
                                      bool refineOnEdgeIntersect,
                                      const string &name,
                                      const bool &debugging=false,
                                      unsigned int tusqhExtraResolveDepth=0,bool Aliasing=false);
        
        virtual bool isItIn(const Polyline &mesh, const list<unsigned int> &faces,
                            const vector<Point3D> &coords) const;

        virtual bool rotateGridMesh(Polyline &input,
									list<RefinementRegion *> &all_reg,
									GeometricTransform &gt);
		
		/*virtual void generateGridFromOctree(const unsigned short &rl, 
                                              Polyline &input,
                                              const string &name);*/
		
        virtual void generateGridMesh(Polyline &input);

        // TUSQH coarse-grid heuristic: when the quadtree starts with
        // very few root cells (typical of large polylines like
        // Chesapeake Bay at level 0), the default TUSQH sample grid
        // (-N 2, only 4 samples per cell) can undersample a cell of
        // size ~50x50 and classify the whole cell as AllInside or
        // AllOutside, producing very few subdivisions. To avoid this
        // we pre-refine the quadtree uniformly to a base level using
        // `generateQuadtreeMesh`, which uses the existing
        // `IntersectionsVisitor` (i.e. legacy edge-intersect refinement)
        // and produces enough small cells for the TUSQH winding
        // criterion to be meaningful. The heuristic is a no-op when
        // the initial grid is already fine enough OR the polyline is
        // simple (few segments, like the unit-square regression cases).
        // Returns true if pre-refinement was applied.
        virtual bool preRefineForTusqh(Polyline &input,
                                       unsigned int maxDepth,
                                       list<RefinementRegion *> &all_reg,
                                       const string &name,
                                       const bool &debugging,
                                       bool Aliasing,
                                       unsigned short coarseThreshold=2,
                                       unsigned short baseLevel=3,
                                       unsigned int minSegmentsForTrigger=100);

        virtual void computeVolumeFractions(Polyline &input, unsigned int sampleSize);

        // TUSQH §3.3 (sub-cell volume fractions). Computes the mean
        // winding number of an s x s sample grid inside a fictitious
        // cell centered on each quadtree vertex and edge. Results are
        // stored in MeshPoint::mSubcellVolumeFraction and in
        // `mEdgeSubcellVF` respectively.
        virtual void computeSubcellVolumeFractions(Polyline &input,
                                                   unsigned int sampleSize,
                                                   double joinThreshold,
                                                   const std::string &name="");

        // TUSQH §3.4 (anti-aliasing - archipelago resolution).
        // 1) compute sub-cell VFs (caller has usually done this),
        // 2) find connected components of the quadtree,
        // 3) for each pair of components connected by an interior edge
        //    (sub-cell VF >= joinThreshold) mark the edge as a "bridge",
        // 4) join the two components by flipping the adjacent quads
        //    along the bridge edge (1-to-5 split pattern), and
        // 5) drop any component that ends up with fewer than
        //    `minComponentCells` leaf quads.
        virtual void resolveArchipelagos(Polyline &input,
                                         unsigned int sampleSize,
                                         double joinThreshold,
                                         unsigned int minComponentCells,
                                         const string& name);

        // TUSQH §3.4 step 3 helper: perform a 1-to-5 split on a quad
        // `q` whose edge at index `bridgeEdgeIdx` (0..3) is the bridge
        // edge. The original quad is REPLACED by 4 sub-quads (from
        // SplitVisitor's standard 1-to-4 split) plus 1 new "bridge"
        // quad that extends OUTWARD from the bridge edge by
        // perpendicular distance H/sampleSize.
        //
        // On success returns a vector of 5 new Quadrants (4 interior
        // sub-quads + 1 bridge quad) with sequential q_ids starting
        // from `nextQIdx`. `nextQIdx` is updated to one past the last
        // assigned q_id.
        //
        // Side effects (the function does these):
        //   - adds new MeshPoints (2 exterior bridge corners + the 4
        //     mid-edges + 1 center from SplitVisitor's 1-to-4 split)
        //     to `points`,
        //   - adds/updates MapEdges entries for the new sub-quads and
        //     the bridge quad,
        //   - adds sub-cell VF placeholder entries for the new edges
        //     so the next `computeSubcellVolumeFractions` call will
        //     find them.
        //
        // On failure returns an empty vector (e.g. degenerate quad
        // where H is too small). The caller must then NOT remove `q`
        // from Quadrants.
        virtual std::vector<Quadrant> bridgeSplitAtEdge(
            Quadrant &q,
            unsigned int bridgeEdgeIdx,
            unsigned int sampleSize,
            unsigned int &nextQIdx);

        // TUSQH §3.4 helper: returns the unit outward direction from the
        // bridge edge of quad `q` (i.e. the direction in which a bridge
        // quad should extend OUTWARD from `edge`).
        //
        // Uses the quad centroid (mean of the 4 corners) as the
        // "interior" reference, so the outward direction is
        // `-normalize(centroid - edge_midpoint)`. This is invariant
        // under rotation and robust to non-convex quads (unlike the
        // midpoint-of-opposite-edge heuristic, which can fail on
        // rotated or non-convex quads).
        //
        // Returns the zero Point3D if `q` is degenerate (centroid
        // coincides with the edge midpoint, |H| < 1e-12).
        virtual Point3D computeExteriorDirection(
            const Quadrant &q,
            unsigned int bridgeEdgeIdx,
            const vector<MeshPoint> &points) const;

        // TUSQH §3.4 helper: decide whether a boundary edge of the
        // quadtree lies on the actual domain boundary (so its
        // "exterior" side is truly outside the geometry, not just
        // outside some quadtree component).
        //
        // Strategy: sample `s x s` points in the rectangular
        // fictitious cell of the edge but restricted to the +dir_ext
        // half-plane. If the mean winding number of those points is
        // STRICTLY below `joinThreshold`, the edge is on the global
        // domain boundary and should NOT be used as a bridge.
        //
        // Returns true ⇔ the edge should be EXCLUDED from the bridge
        // candidate set (i.e. it sits on the global domain boundary).
        virtual bool isEdgeOnDomainBoundary(
            const QuadEdge &edge,
            const Quadrant &adjacentQuad,
            const Polyline &input,
            unsigned int sampleSize,
            double joinThreshold) const;

        // Paper-faithful helper: returns true iff a Quadrant is an
        // "interior" cell of the cubical complex, i.e. it has positive
        // volume fraction (AllInside) or is on the boundary (Mixed).
        // AllOutside cells return false and are skipped by:
        //   - the bridge-joining BFS (they don't count as part of
        //     any component but they DO remain in `Quadrants` so that
        //     MapEdges retains the edges to neighbouring cells),
        //   - `saveOutputMesh` (excluded from the final mesh output),
        //   - any other downstream consumer that needs the
        //     geometric mesh (linkElementsToNodes,
        //     detectFeatureQuadrants, detectInsideNodes, etc.).
        static bool isInteriorCell(const Quadrant &q);

        virtual void detectFeatureQuadrants(Polyline &input);
		
        virtual void linkElementsToNodes();

        virtual void computeNodeMaxDist();

        virtual void detectInsideNodes(Polyline &input);

		virtual void removeOnSurface(Polyline &input);
        
        virtual void removeOnSurfaceSafe(Polyline &input);
		
        virtual void applySurfacePatterns(Polyline &input);

        virtual void shrinkToBoundary(Polyline &input);

        virtual unsigned int saveOutputMesh(const shared_ptr<FEMesh> &mesh, bool decoration=false);
		
        virtual unsigned int saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                            vector<MeshPoint> &points,
                                            list<Quadrant> &elements,
                                            const bool &debugging=false,
                                            const list<Point3D> &extra_pts=list<Point3D> ());
        
        virtual unsigned int saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                            vector<MeshPoint> &points,
                                            vector<Quadrant> &elements,
                                            const bool &debugging=false,
                                            const list<Point3D> &extra_pts=list<Point3D> ());


        virtual void projectCloseToBoundaryNodes(Polyline &input);


		
	protected:
		
		vector<MeshPoint> points;
		vector<Quadrant> Quadrants;
        //Map that for each edge saves its mid point index (0 by default).
        map<QuadEdge, EdgeInfo> MapEdges;
		list<RefinementRegion *> regions;

        unsigned int mSampleSize;

        // Per-edge sub-cell volume fraction. Kept in a side map to
        // avoid disturbing the on-disk format of EdgeInfo (used by
        // Services::WriteQuadtreeMesh / ReadQuadMesh). Keys are the
        // canonical (sorted) edge endpoints, values are the s*s grid
        // subgrid VF computed by computeSubcellVolumeFractions.
        std::map<QuadEdge, EdgeSubcellVFData> mEdgeSubcellVF;



	};
}
#endif
