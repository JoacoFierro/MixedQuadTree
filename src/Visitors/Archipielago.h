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
#include "../EdgeInfo.h"

#include <string>

#ifndef Archipielago_h
#define Archipielago_h 1

using Clobscode::MeshPoint;
using Clobscode::QuadEdge;
using Clobscode::Point3D;
using Clobscode::Quadrant;
using std::vector;
using std::list;
using std::set;
using std::pair;

struct TemplateInfo{
    QuadEdge edge;
    const Quadrant* quad1 = nullptr;
    const Quadrant* quad2 = nullptr;
};

struct FixWithQuad{
    unsigned int sharedVertex;
    TemplateInfo qTemplate;
    const Quadrant *qQuadrant;
};

struct FixWithTemp{
    unsigned int sharedVertex;
    TemplateInfo Template1;
    TemplateInfo Template2;

};
namespace Clobscode
{
    class Archipielago {
        public:
        
            virtual void setQuadrant(list<Quadrant> &Quadrants, list<Quadrant> BackQuadrants);

            virtual void setPoints (vector<MeshPoint> &Meshpoints, vector<MeshPoint> BackPoints);

            virtual void setActualIndex(unsigned int &idx);

            virtual void setMapEdges(map<QuadEdge, EdgeInfo> MapEdges, map<QuadEdge, EdgeInfo> BackMapEdges);

            virtual void setInput(Polyline &ply);

            virtual void getTemplatesToFix();

            virtual void createTemplates(const QuadEdge edge, unsigned int maxDepth);

            virtual void createFixTempQuad(const QuadEdge edge, unsigned int maxDepth,const Quadrant q);

            virtual void createFixWithTemp(unsigned int sharedPoint,
                              unsigned int p1a,
                              unsigned int p1b,
                              unsigned int p2a,
                              unsigned int p2b,
                              unsigned int maxDepth);

            virtual bool SubSampling(const QuadEdge edge);

            virtual void getOutsideEdges( unsigned int maxDepth);

            virtual void fixTemplates(unsigned int maxDepth);

            virtual QuadEdge SelectEdge(QuadEdge Tedge,QuadEdge Qedge1,QuadEdge Qedge2,unsigned int shared);

            virtual unsigned int getIndexWithSharedVertex(Quadrant q,unsigned int shared,QuadEdge edge);

            virtual bool pointInsideTemplate(const Quadrant &q,const Point3D &P) const;

            virtual unsigned int getTemplatesInsertados() const { return TemplatesInsertados; }

            virtual void saveFixInformation(const std::string& filename) const;

        protected:
            list<Quadrant> *Quadrants;
            vector<MeshPoint> *points;
            map<QuadEdge, EdgeInfo> MapEdges;

            list<Quadrant> BackQuadrants;
            vector<MeshPoint> BackPoints;
            map<QuadEdge, EdgeInfo> BackMapEdges;

            Polyline *ply;
            unsigned int *CurrentQuadIndex;

            std::map<unsigned int, std::vector<const Quadrant*>> mapVertexQuadrants;
            vector<FixWithQuad> ToFixWithQuads;
            vector<FixWithTemp> ToFixWithTems;
            vector<TemplateInfo> TemplesAdded;

            unsigned int TemplatesInsertados = 0;

    };
}

#endif
