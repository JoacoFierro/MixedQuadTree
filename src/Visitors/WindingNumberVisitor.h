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
 * @file WindingNumberVisitor.h
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief Visitor for computing winding numbers and volume fractions on QuadTree
 **/

#ifndef WindingNumberVisitor_h
#define WindingNumberVisitor_h

#include "Visitor.h"
#include <vector>

using std::vector;

namespace Clobscode
{
    class Polyline;
    class MeshPoint;
    class Quadrant;

    class WindingNumberVisitor : public Visitor
    {
    public:

        WindingNumberVisitor(unsigned int s);

        virtual ~WindingNumberVisitor();

        virtual bool visit(Quadrant *q) override;

        void setPolyline(Polyline *geo);

        void setPoints(vector<MeshPoint> *mp);

    private:

        void computePostOrder(Quadrant *q);

        unsigned int mSampleSize;

        Polyline* mPolyline;

        vector<MeshPoint>* mPoints;
    };
}

#endif