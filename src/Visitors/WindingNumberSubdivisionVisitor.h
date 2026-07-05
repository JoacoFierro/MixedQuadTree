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
  * @file WindingNumberSubdivisionVisitor.h
  * @author Claudio Lobos, Fabrice Jaillet
  * @version 0.1
  * @brief Visitor that decides whether a quadrant must be subdivided based
  *        on the TUSQH criterion: a cell is subdivided only when its s x s
  *        sample points do not all agree on being inside (winding number > 0)
  *        or outside (winding number == 0) the input domain.
  *
  *        Reference:
  *          "Tusqh: Topological Control of Volume-Fraction Meshes Near Small
  *           Features and Dirty Geometry" (Bracci et al.).
  **/
#ifndef WindingNumberSubdivisionVisitor_h
#define WindingNumberSubdivisionVisitor_h 1

#include "Visitor.h"
#include <vector>

using std::vector;

namespace Clobscode
{
    class Polyline;
    class MeshPoint;

    class WindingNumberSubdivisionVisitor : public Visitor
    {
    public:

        // s : number of samples per axis (s x s grid of sample points).
        // refineOnEdgeIntersect : when true, also force refinement of cells
        //     that intersect any Polyline edge regardless of the winding
        //     ambiguity (this preserves the legacy behaviour for users
        //     that want pure edge-based refinement). Defaults to false.
        WindingNumberSubdivisionVisitor(unsigned int s,
                                        bool refineOnEdgeIntersect = false);

        virtual ~WindingNumberSubdivisionVisitor();

        virtual bool visit(Quadrant *q) override;

        void setPolyline(Polyline *geo);

        void setPoints(vector<MeshPoint> *mp);

        // Optional: if provided, a quadrant that intersects the Polyline
        // geometrically will always be refined (legacy criterion).
        // The IntersectionsVisitor::select_edges=true is what is used here.
        // Pass nullptr to disable.
        void setIntersectionsVisitor(class IntersectionsVisitor *iv);

    private:

        unsigned int mSampleSize;
        bool mRefineOnEdgeIntersect;
        Polyline *mPolyline;
        vector<MeshPoint> *mPoints;
        IntersectionsVisitor *mIntersectionsVisitor;
    };
}

#endif
