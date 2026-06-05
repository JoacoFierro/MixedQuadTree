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
  * @file WindingNumberVisitor.cpp
  * @author Claudio Lobos, Fabrice Jaillet
  * @version 0.1
  * @brief Visitor for computing winding numbers and volume fractions on QuadTree
  **/

#include "WindingNumberVisitor.h"
#include "../Quadrant.h"
#include "../Polyline.h"

namespace Clobscode
{

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    WindingNumberVisitor::WindingNumberVisitor(unsigned int s)
        :mSampleSize(s), mPolyline(nullptr), mPoints(nullptr), mQuadrants(nullptr)
    {
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    WindingNumberVisitor::~WindingNumberVisitor()
    {
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::setPolyline(Polyline *geo)
    {
        mPolyline = geo;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::setPoints(vector<MeshPoint> *mp)
    {
        mPoints = mp;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::setQuadrants(vector<Quadrant> *quadrants)
    {
        mQuadrants = quadrants;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool WindingNumberVisitor::visit(Quadrant *q)
    {
        q->setSampleSize(mSampleSize);

        if (mQuadrants == nullptr) {
            std::cerr << "Error: WindingNumberVisitor::visit() - mQuadrants is null\n";
            return false;
        }

        computePostOrder(q);
        return true;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::computePostOrder(Quadrant *q)
    {
        const auto& subs = q->getSubElements();

        if (subs.size() > 1) {
            q->mNeedsInheritance = true;
            return;
        }

        unsigned int s = mSampleSize;
        vector<double> wn(s * s);

        for (unsigned int i = 0; i < s; i++) {
            for (unsigned int j = 0; j < s; j++) {
                Point3D sample = q->getSamplePoint(i, j, *mPoints);
                wn[i * s + j] = mPolyline->windingNumber(sample);
            }
        }

        q->computeVolumeFraction(wn);
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::computeInheritance(Quadrant *q)
    {
        const auto& subs = q->getSubElements();

        if (subs.size() <= 1) {
            return;
        }

        if (!q->mNeedsInheritance) {
            return;
        }

        double sum_vf = 0.0;
        int count_children_with_vf = 0;

        for (const auto& child_sub : subs) {
            for (auto& child : *mQuadrants) {
                const vector<unsigned int>& child_points = child.getPointIndex();
                if (child_points.size() == 4 && child_sub.size() == 4) {
                    bool match = (child_points[0] == child_sub[0] &&
                                  child_points[1] == child_sub[1] &&
                                  child_points[2] == child_sub[2] &&
                                  child_points[3] == child_sub[3]);
                    if (match) {
                        if (child.hasVolumeFraction()) {
                            sum_vf += child.getVolumeFraction();
                            count_children_with_vf++;
                        }
                        computeInheritance(&child);
                        break;
                    }
                }
            }
        }

        if (count_children_with_vf > 0) {
            q->mVolumeFraction = sum_vf / count_children_with_vf;
            q->mHasVolumeFraction = true;
        }
        q->mNeedsInheritance = false;
    }

}