/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

 Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

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
* @file MeshPoint.cpp
* @author Claudio Lobos, Fabrice Jaillet
* @version 0.1
* @brief
**/

#include "MeshPoint.h"

namespace Clobscode
{

    //MeshPoint::MeshPoint():outsidechecked(false), projected(false),feature(false),inside(false),maxdistance(std::numeric_limits<double>::max()){
    MeshPoint::MeshPoint(): state(STATEMASK),maxdistance(std::numeric_limits<double>::max()),
                            mSubcellSampleSize(0), mSubcellVolumeFraction(0.0),
                            mHasSubcellVolumeFraction(false), mSubcellIsInterior(false){
        //we assume that every point is outside by default.(inside)
        //checking if a point is outside or not is a very expensive
        //operation, so we try to do it only once (outsidechecked)
	}

    //MeshPoint::MeshPoint(const Point3D &p):outsidechecked(false), projected(false),feature(false),inside(false),maxdistance(std::numeric_limits<double>::max()){
    MeshPoint::MeshPoint(const Point3D &p):point(p),state(STATEMASK),maxdistance(std::numeric_limits<double>::max()),
                                  mSubcellSampleSize(0), mSubcellVolumeFraction(0.0),
                                  mHasSubcellVolumeFraction(false), mSubcellIsInterior(false){
//          point = p;
    }

	MeshPoint::~MeshPoint(){

	}

    void MeshPoint::computeSubcellVolumeFraction(const std::vector<double>& wn) {
        mSubcellWindingNumbers = wn;
        double sum = 0.0;
        for (double v : wn) sum += v;
        mSubcellVolumeFraction = wn.empty() ? 0.0 : sum / static_cast<double>(wn.size());
        mHasSubcellVolumeFraction = true;
    }

    std::ostream& operator<<(std::ostream& o, const MeshPoint &p){
		o << p.getPoint();
		return o;
	}
	
}
