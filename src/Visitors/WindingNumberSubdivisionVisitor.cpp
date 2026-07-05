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
  * @file WindingNumberSubdivisionVisitor.cpp
  * @author Claudio Lobos, Fabrice Jaillet
  * @version 0.1
  * @brief TUSQH-style subdivision criterion (see WindingNumberSubdivisionVisitor.h).
  *
  * Subdivision rule (TUSQH):
  *   - Compute winding number wn(p) for each of the s x s sample points of
  *     the cell.
  *   - If all samples are strictly positive (wn>0) -> cell is "AllInside",
  *     no subdivision.
  *   - If all samples are zero (wn==0) -> cell is "AllOutside", no
  *     subdivision.
  *   - If samples are mixed -> cell is "Mixed", must be subdivided.
  *   - Optional legacy behaviour: if a cell intersects any input edge
  *     geometrically, force subdivision regardless of the winding
  *     criterion (this allows reproducing the original edge-based
  *     refinement for benchmarking).
  *
  * The visitor also stores the s x s winding numbers and the volume
  * fraction on the quadrant, so no extra pass is required.
  **/

#include "WindingNumberSubdivisionVisitor.h"
#include "../Quadrant.h"
#include "../Polyline.h"
#include "IntersectionsVisitor.h"

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    WindingNumberSubdivisionVisitor::WindingNumberSubdivisionVisitor(unsigned int s,
                                                                     bool refineOnEdgeIntersect)
        :mSampleSize(s),
         mRefineOnEdgeIntersect(refineOnEdgeIntersect),
         mPolyline(nullptr),
         mPoints(nullptr),
         mIntersectionsVisitor(nullptr)
    {
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    WindingNumberSubdivisionVisitor::~WindingNumberSubdivisionVisitor()
    {
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberSubdivisionVisitor::setPolyline(Polyline *geo)
    {
        mPolyline = geo;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberSubdivisionVisitor::setPoints(vector<MeshPoint> *mp)
    {
        mPoints = mp;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberSubdivisionVisitor::setIntersectionsVisitor(IntersectionsVisitor *iv)
    {
        mIntersectionsVisitor = iv;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool WindingNumberSubdivisionVisitor::visit(Quadrant *q)
    {
        if (mPolyline == nullptr || mPoints == nullptr) {
            return false;
        }

        q->setSampleSize(mSampleSize);

        const unsigned int s = mSampleSize;

        // Legacy criterion (optional): if the cell geometrically intersects
        // any input edge, force subdivision. By default this is OFF so that
        // we strictly follow the TUSQH criterion.
        if (mRefineOnEdgeIntersect && mIntersectionsVisitor != nullptr) {
            // Use the IntersectionsVisitor with the parent's intersected
            // edges as candidate list. Empty list means "check all edges"
            // but that is expensive, so we only re-run when we already
            // have intersected edges to test.
            if (!q->getIntersectedEdges().empty()) {
                // We have a smaller candidate list, run the visitor on it.
                std::list<unsigned int> local(q->getIntersectedEdges());
                std::vector<Point3D> coords(2);
                coords[0] = (*mPoints)[q->getPointIndex()[0]].getPoint();
                coords[1] = (*mPoints)[q->getPointIndex()[2]].getPoint();
                mIntersectionsVisitor->setPolyline(*mPolyline);
                mIntersectionsVisitor->setEdges(local);
                mIntersectionsVisitor->setCoords(coords);
                // We just need the yes/no answer; the visitor overwrites
                // the quadrant's intersected edges, so do a temporary copy
                // and restore.
                std::list<unsigned int> backup = q->getIntersectedEdges();
                bool hits = q->accept(mIntersectionsVisitor);
                q->setIntersectedEdges(backup);
                if (hits) {
                    q->setWindingState(WindingState::Mixed);
                    // Note: we don't recompute the winding numbers /
                    // volume fraction here. The caller (Mesher) will do
                    // that later on the children once the loop stops.
                    return true;
                }
            }
        }

        // TUSQH main criterion: classify the cell from the s x s
        // winding-number samples.
        std::vector<double> wn_values;
        wn_values.reserve(s * s);

        bool anyPositive = false;
        bool anyZero = false;

        for (unsigned int i = 0; i < s; ++i) {
            for (unsigned int j = 0; j < s; ++j) {
                Point3D sample = q->getSamplePoint(i, j, *mPoints);
                int wn = mPolyline->windingNumber(sample);
                if (wn > 0) {
                    anyPositive = true;
                    wn_values.push_back(1.0);
                } else {
                    anyZero = true;
                    wn_values.push_back(0.0);
                }
            }
        }

        if (!anyPositive && !anyZero) {
            // Shouldn't happen with s>0.
            q->setWindingState(WindingState::Unknown);
            return false;
        }

        if (anyPositive && !anyZero) {
            q->setWindingState(WindingState::AllInside);
            // Reuse the same machinery that WindingNumberVisitor uses so
            // that volume fraction / sample data are available downstream.
            q->computeVolumeFraction(wn_values);
            return false;
        }
        if (!anyPositive && anyZero) {
            q->setWindingState(WindingState::AllOutside);
            q->computeVolumeFraction(wn_values);
            return false;
        }

        // Mixed: must be subdivided per TUSQH. We still record the
        // current wn values so that the user can inspect the heatmap of
        // the cell as it stands now.
        q->setWindingState(WindingState::Mixed);
        q->computeVolumeFraction(wn_values);
        return true;
    }
}
