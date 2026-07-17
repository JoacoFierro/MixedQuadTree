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
* @file PinchVisitor.cpp
* @brief Trivial implementation of PinchVisitor (setters only).
*        Concrete subclasses implement `visit`.
**/

#include "PinchVisitor.h"

namespace Clobscode
{

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    PinchVisitor::PinchVisitor()
        : mPoints(nullptr),
          mQuadrants(nullptr),
          mMapEdges(nullptr),
          mJoinThreshold(0.5)
    {
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    PinchVisitor::~PinchVisitor()
    {
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    void PinchVisitor::setPoints(const std::vector<MeshPoint> *mp)
    {
        mPoints = mp;
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    void PinchVisitor::setQuadrants(const std::vector<Quadrant> *q)
    {
        mQuadrants = q;
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    void PinchVisitor::setMapEdges(const std::map<QuadEdge, EdgeInfo> *me)
    {
        mMapEdges = me;
    }

    //--------------------------------------------------------------------
    //--------------------------------------------------------------------
    void PinchVisitor::setJoinThreshold(double t)
    {
        mJoinThreshold = t;
    }

} // namespace Clobscode
