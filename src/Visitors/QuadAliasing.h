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
* @file QuadAliasing.h
* @author Joaquin Fierro, Felipe
* @version 0.1
* @brief
**/

#ifndef QuadAliasing_h
#define QuadAliasing_h 1

#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <cmath>
#include <math.h>
#include "../Quadrant.h"
#include "../MeshPoint.h"
#include "../QuadEdge.h"
#include "../Polyline.h"
#include <string>


using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::Point3D;
using Clobscode::Quadrant;
using std::vector;
using std::list;
using std::set;
using std::pair;

struct VertexAlias{

    unsigned int sharedVertex;
    const Quadrant *q1;
    const Quadrant *q2;

    std::array<QuadEdge,2> q1_edges;
    std::array<QuadEdge,2> q2_edges;
};

namespace Clobscode
{
    class QuadAliasing {
        public:

            virtual void setQuadrant(list<Quadrant> &Quadrants);

            virtual void setPoints (vector<MeshPoint> &Meshpoints);

            virtual void setActualIndex(unsigned int idx);

            virtual void setInput(const Polyline &ply);

            virtual void printQuadrants();

            virtual void printPoints();

            virtual void printPinches();

            virtual void savePinches(const std::string &filename);

            virtual void getPinches();

            virtual void getAdjacentEdges(const Quadrant &q, unsigned int vertex, QuadEdge &e1, QuadEdge &e2);

            virtual bool edgeHasNeighbour(const QuadEdge &e);

            virtual void QuadrantVertexMap();

            virtual void CreateTemplates();

            virtual void createTemplate(const VertexAlias &pinch);

            virtual unsigned int oppositeVertex(const QuadEdge &e, unsigned int shared);

            virtual bool ComparePoints (const Point3D &p1, const Point3D &p2);

            bool validAngularConfiguration(unsigned int shared,unsigned int A,unsigned int B, unsigned int C,unsigned int D);

        protected:
            list<Quadrant> *Quadrants;
            vector<MeshPoint> *points;
            const Polyline *ply;
            std::map<QuadEdge, std::vector<const Quadrant*>> mapQuadrant;
            std::map<unsigned int, std::vector<const Quadrant*>> mapVertexQuadrants;
            vector<VertexAlias> Pinches;

            unsigned int CurrentQuadIndex;
    };
}

#endif
