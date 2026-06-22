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
  *
  * Based on Tusqh paper: "The volume fraction of each maximal-dimension cell
  * is computed as the average of the winding numbers of its sample points."
  **/

#include "WindingNumberVisitor.h"
#include "../Quadrant.h"
#include "../Polyline.h"

namespace Clobscode
{

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    WindingNumberVisitor::WindingNumberVisitor(unsigned int s)
        :mSampleSize(s), mPolyline(nullptr), mPoints(nullptr)
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
    bool WindingNumberVisitor::visit(Quadrant *q)
    {
        q->setSampleSize(mSampleSize);
        computePostOrder(q);
        return true;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void WindingNumberVisitor::computePostOrder(Quadrant *q)
    {
        vector<double> wn_values;
        unsigned int s = mSampleSize;

        for (unsigned int i = 0; i < s; ++i) {
            for (unsigned int j = 0; j < s; ++j) {
                Point3D sample = q->getSamplePoint(i, j, *mPoints);
                double wn = mPolyline->windingNumber(sample);
                wn_values.push_back(wn);
            }
        }

        if (!wn_values.empty()) {
            q->computeVolumeFraction(wn_values);
        }
    }

}