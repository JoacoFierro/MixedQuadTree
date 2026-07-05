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
                                  unsigned int tusqhExtraResolveDepth=0);

        virtual std::shared_ptr<FEMesh> refineMesh(Polyline &input, const unsigned short &rl,
                                  const string &name, list<unsigned int> &roctli,
                                  list<RefinementRegion *> &all_reg,
                                  GeometricTransform &gt, const unsigned short &minrl,
                                  const unsigned short &omaxrl, const bool &debugging,
                                  unsigned int sampleSize=2, bool decoration=false,
                                  bool useTusqh=false,
                                  unsigned int tusqhSampleSize=2,
                                  bool refineOnEdgeIntersect=false,
                                  unsigned int tusqhExtraResolveDepth=0);

        
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
                                          unsigned int new_q_idx=0);

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
                                      unsigned int tusqhExtraResolveDepth=0);
        
        virtual bool isItIn(const Polyline &mesh, const list<unsigned int> &faces,
                            const vector<Point3D> &coords) const;

        virtual bool rotateGridMesh(Polyline &input,
									list<RefinementRegion *> &all_reg,
									GeometricTransform &gt);
		
		/*virtual void generateGridFromOctree(const unsigned short &rl, 
                                              Polyline &input,
                                              const string &name);*/
		
        virtual void generateGridMesh(Polyline &input);

        virtual void computeVolumeFractions(Polyline &input, unsigned int sampleSize);

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



	};
}
#endif
