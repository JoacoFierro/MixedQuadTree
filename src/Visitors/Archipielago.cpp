/*
 <Mix-mesher: region type. This program generates a mixed-elements mesh>
 
 Copyright (C) <2013,2017>  <Claudio Lobos>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/gpl.txt>
 */
#include <array>
#include <map>
#include "Archipielago.h"
#include "../GridMesher.h"
#include "QuadAliasing.h"

 namespace Clobscode {

    void Archipielago::setQuadrant( list<Quadrant> &quadrants,list <Quadrant> BackQuadrants){
        this->Quadrants = &quadrants;
        this->BackQuadrants = BackQuadrants;
    }

    void Archipielago::setPoints( vector<MeshPoint> &Meshpoints, vector<MeshPoint> BackPoints){
        this->points = &Meshpoints;
        this->BackPoints = BackPoints;
    }

    void Archipielago::setMapEdges(map<QuadEdge, EdgeInfo> MapEdges, map<QuadEdge, EdgeInfo> BackMapEdges){
        this->MapEdges = MapEdges;
        this->BackMapEdges = BackMapEdges;
    }

    void Archipielago::setInput(Polyline &ply){
        this->ply = &ply;
    }

    void Archipielago::setActualIndex(unsigned int &idx){
        this->CurrentQuadIndex = &idx;
    }
    
    void Archipielago::getTemplatesToFix(){
        bool isTemplate;
        unsigned int sharedVerted, sharedVertex , id0, id1 , T1_opposite, T2_opossite;
        cout<<"Cantidad de templates archipielagos obtenidos: " << TemplesAdded.size() <<endl;
        for (size_t i = 0; i < TemplesAdded.size(); ++i) {
            const TemplateInfo &T1 = TemplesAdded[i];
            id0 = T1.edge[0];
            id1 = T1.edge[1];
            for(const Quadrant &q: *Quadrants){
                unsigned int count = 0; 
                isTemplate = q.getIsTemplate();
                const vector<unsigned int> idx = q.getPointIndex();
                for (unsigned int i=0; i!=idx.size();i++){
                    if (id0 == idx[i] || id1 == idx[i]){
                        sharedVerted = idx[i];
                        count++;
                    }
                }
                if (count==1){
                    isTemplate = q.getIsTemplate();
                    if(isTemplate == false){
                        FixWithQuad Fa;
                        Fa.sharedVertex = sharedVerted;
                        Fa.qTemplate = T1;
                        Fa.qQuadrant = &q;
                        ToFixWithQuads.push_back(Fa);
                    }
                }
            }
            for (size_t j = i + 1; j < TemplesAdded.size(); ++j) {
                const TemplateInfo &T2 = TemplesAdded[j];
                if (T2.edge[0] == id0 || T2.edge[0] == id1)
                    sharedVertex = T2.edge[0];
                else if (T2.edge[1] == id0 || T2.edge[1] == id1)
                    sharedVertex = T2.edge[1];
                else
                    continue;

                T1_opposite = sharedVertex ==  id0 ? id1: id0;
                T2_opossite = sharedVertex == T2.edge[0]? T2.edge[1] : T2.edge[0];
                
                Point3D PA = (*points)[sharedVertex].getPoint();
                Point3D PB = (*points)[T1_opposite].getPoint();
                Point3D PC = (*points)[T2_opossite].getPoint();

                Point3D v1 = PB - PA;
                Point3D v2 = PC - PA;

                v1.normalize();
                v2.normalize();

                double dot = v1 * v2; 
                const double eps = 1e-6;

                if (abs(dot) < eps) continue;

                FixWithTemp ft;
                ft.sharedVertex = sharedVertex;
                ft.Template1 = T1;
                ft.Template2 = T2;

                ToFixWithTems.push_back(ft);
            }
        }
        cout<<"ToFixWithQuads: "<< ToFixWithQuads.size() << endl;
        cout<<"ToFixWithTems: " << ToFixWithTems.size() <<endl;
    }

    void Archipielago::createFixWithTemp(unsigned int sharedPoint,
                                         unsigned int Id1a,
                                         unsigned int Id1b,
                                         unsigned int Id2a,
                                         unsigned int Id2b,
                                         unsigned int maxDepth){

        // Obtenemos los puntos 

        Point3D S   = (*points)[sharedPoint].getPoint();
        Point3D p1a = (*points)[Id1a].getPoint();
        Point3D p1b = (*points)[Id1b].getPoint();
        Point3D p2a = (*points)[Id2a].getPoint();
        Point3D p2b = (*points)[Id2b].getPoint();
        
        double eps = 1e-6, width = 0.7;
        bool confAA = true; 
        bool isVertical = abs(p1a.Y() - p1b.Y()) < eps;
               
        Point3D dir;
        if(isVertical) // Es vertical la configuracion
            if (abs(p1a.X() - p2a.X()) < eps){
                dir = p1a -  p2a;
            } else{
                dir = p1a -  p2b;
                confAA = false;
        } else {                        // Es Horizontal  la configuracion
            if (abs(p1a.Y() - p2a.Y()) < eps) { 
            } else {
                dir = p1a -  p2b;
                confAA = false;
            }
        }

        dir.normalize();
        Point3D n (-dir.Y(), dir.X(),0);
        n.normalize();

        Point3D I = S + n * width;
        Point3D D = S - n * width;

        // Verificar si esta dentrode un template o cuadrante
        bool insertPointI = true, insertPointD = true;
        unsigned int idI,idD; 
        
        QuadEdge EdgeVerify;
        for(TemplateInfo &ToRepair: TemplesAdded){
            EdgeVerify = ToRepair.edge;
            Point3D A = (*points)[EdgeVerify[0]].getPoint();
            Point3D B = (*points)[EdgeVerify[1]].getPoint();
            if(isVertical){
                if( min(A.X(),B.X()) < I.X() && I.X() < max(A.X(),B.X()))insertPointI = false;
                if( min(A.X(),B.X()) < D.X() && D.X() < max(A.X(),B.X()))insertPointD = false;  
            } else{
                if( min(A.Y(),B.Y()) < I.Y() && I.Y() < max(A.Y(),B.Y()))insertPointI = false;
                if( min(A.Y(),B.Y()) < D.Y() && D.Y() < max(A.Y(),B.Y()))insertPointD = false;  
            }
        }

        if(insertPointI){
            idI=points->size();
            points->push_back(MeshPoint(I));
        }
        if(insertPointD){
            idD=points->size();
            points->push_back(MeshPoint(D));
        }

        //Determinanmos la correspondencia entre puntos

        double da = I.distance(p1a);
        double db = I.distance(p1b);

        if(da < db ? insertPointI :insertPointD ){
            vector<unsigned int> QuadA = { (da < db ? idI : idD) ,  Id1a, sharedPoint, (confAA ? Id2a: Id2b)}; 
            Quadrant q1 (QuadA,maxDepth, (*CurrentQuadIndex)++);
            q1.setIsTemplate(true);
            Quadrants->push_back(q1);
            TemplatesInsertados++;
        }
        if(da < db ? insertPointD :insertPointI ){
            vector<unsigned int> QuadB = { (da < db ?idD : idI),  (confAA ? Id2b: Id2a), sharedPoint, Id1b }; 
            Quadrant q2 (QuadB,maxDepth, (*CurrentQuadIndex)++);
            q2.setIsTemplate(true);
            Quadrants->push_back(q2);
            TemplatesInsertados++;
        }
    }


    void Archipielago::createFixTempQuad(const QuadEdge edge, unsigned int maxDepth,const Quadrant q){
        unsigned int idp0 = edge[0];
        unsigned int idp1 = edge[1];

        const Point3D p0 = (*points)[idp0].getPoint();
        const Point3D p1 = (*points)[idp1].getPoint();

        Point3D dir = p1 - p0;
        double L = dir.Norm();
        double width = 0.45 * L;

        dir.normalize();
        Point3D n(-dir.Y(), dir.X(), 0.0);
        n.normalize();

        Point3D mid = (p0 + p1) * 0.5;
        Point3D test = mid + n * (0.1 * (p1-p0).Norm());
        if (q.pointInside(*points, test)) { // Si la normal esta adentro del cuadrante
            n = -n;
        }

        Point3D bevel0a = dir + n;
        Point3D bevel1a = -dir + n;

        bevel0a.normalize();
        bevel1a.normalize();

        Point3D A = p0 + bevel0a*width;
        Point3D C = p1 + bevel1a*width;

        bool insertPointA  = true;
        bool insertPointC  = true;
        unsigned int idA,idC;

        QuadAliasing qa;
        for (unsigned int i = 0; i < points->size(); ++i) {
            const Point3D &p = (*points)[i].getPoint();
            if (qa.ComparePoints(A,p)){
                idA = i;
                insertPointA = false;
            }
            if (qa.ComparePoints(C,p)){
                idC = i;
                insertPointC = false;
            }
            if (!insertPointA && !insertPointC)
                break;
        }
        bool insertT = true;
        Point3D midt = (A + C) * 0.5;

        double delta = (p1 - p0).Norm() * 0.05; 
        Point3D exterior = (p0 + p1) * 0.5 + n * delta;

        for (const Quadrant &q: *Quadrants){
            if(q.getIsTemplate()){

                if(pointInsideTemplate(q,exterior))insertT= false;
            }
        }

        for (auto templateEdge:  TemplesAdded){
            const Quadrant* quad1 = templateEdge.quad1;
            const Quadrant* quad2 = templateEdge.quad2;
            QuadEdge e = templateEdge.edge;
            Point3D p1 = (*points)[e[0]].getPoint();
            Point3D p2 = (*points)[e[1]].getPoint();

            // if (quad2->pointInside(*points, midt) || quad2->pointInside(*points, midt)) insertT = false;
            
            if (abs(p1.X() - p2.X()) < 1E-8){
                if(abs(p1.X() - midt.X()) < 1E-8){
                    if (midt.Y() > min(p1.Y(),p2.Y()) && midt.Y() < max(p1.Y(),p2.Y())) insertT = false;
                }
            } else{
                if(abs(p1.Y() - midt.Y()) < 1E-8){
                    if (midt.X() > min(p1.X(),p2.X()) && midt.X() < max(p1.X(),p2.X()))insertT = false;
                }
            }
        }
        
        if(insertT){

            if(insertPointA){
                idA=points->size();
                points->push_back(MeshPoint(A));
            }
            if(insertPointC){
                idC=points->size();
                points->push_back(MeshPoint(C));
            }
            
            vector<unsigned int> QuadAC = {idp0,idA,idC,idp1};
            Quadrant q1 (QuadAC,maxDepth, (*CurrentQuadIndex)++);
            q1.setIsTemplate(true);
            QuadEdge e1(idp0,idp1);
            q1.setTemplateEdge(e1);
            Quadrants->push_back(q1);
            TemplatesInsertados++;
        }

    }

    void Archipielago::createTemplates(const QuadEdge edge, unsigned int maxDepth){

        const Point3D p0  = (BackPoints)[edge[0]].getPoint();
        const Point3D p1  = (BackPoints)[edge[1]].getPoint();

        unsigned int idp0;
        unsigned int idp1;

        bool insertPointsP0 = true;
        bool insertPointsP1 = true;

        QuadAliasing qa;
        for (unsigned int i = 0; i < points->size(); i++) {
            const Point3D &p = (*points)[i].getPoint();
            if (qa.ComparePoints(p0,p)){
                idp0 = i;
                insertPointsP0 = false;
            }
            if (qa.ComparePoints(p1,p)){
                idp1 = i;
                insertPointsP1 = false;
            }
            if (!insertPointsP0 && !insertPointsP1)
                break;
        }

        if(insertPointsP0){
            idp0=points->size();
            points->push_back(MeshPoint(p0));
        }
        if(insertPointsP1){
            idp1=points->size();
            points->push_back(MeshPoint(p1));
        }

        Point3D dir = p1 - p0;
        double L = dir.Norm();
        dir.normalize();

        Point3D n(-dir.Y(), dir.X(), 0.0);
        n.normalize();

        double width = 0.45 * L;

        Point3D bevel0a, bevel0b, bevel1a, bevel1b;
        Point3D A,B,C,D;

        bevel0a = dir + n; //A
        bevel1a = -dir + n; // B

        bevel0b = dir - n; // C 
        bevel1b = -dir - n; // D

        bevel0a.normalize();
        bevel1a.normalize();

        bevel0b.normalize();
        bevel1b.normalize();

        A = p0 + bevel0a*width;
        C = p1 + bevel1a*width;

        B = p0 + bevel0b*width;
        D = p1 + bevel1b*width;

        bool InsertA = true, InsertB = true, InsertC = true, InsertD = true;

        for (const Quadrant &q : *Quadrants) {
            // if (q.getIsTemplate()){
                if (q.pointInside(*points,A)){
                    InsertA = false;
                }
                if (q.pointInside(*points,A)) {
                    InsertB = false;
                }
                if (q.pointInside(*points,C)) {
                    InsertC = false;
                }
                if (q.pointInside(*points,D)) {
                    InsertD = false;
                }
            // }
        }

        if ((InsertA && InsertC) && (InsertB && InsertD)){
            TemplateInfo ti;
            QuadEdge e1(idp0,idp1);
            ti.edge = e1;
            //Fist Quadrant
            unsigned int idA=points->size();
            points->push_back(MeshPoint(A));
            unsigned int idC=points->size();
            points->push_back(MeshPoint(C));

            vector<unsigned int> QuadAC = {idp0,idA,idC,idp1};
            Quadrant q1 (QuadAC,maxDepth, (*CurrentQuadIndex)++);
            q1.setIsTemplate(true);
            
            Quadrants->push_back(q1);
            Quadrant *ptrQ1 = &Quadrants->back();
            ti.quad1 = ptrQ1;
            TemplatesInsertados++;
  

            //Second Quadrant
            unsigned int idB=points->size();
            points->push_back(MeshPoint(B));
            unsigned int idD=points->size();
            points->push_back(MeshPoint(D));
            
            vector<unsigned int> QuadBD = {idp0,idB,idD,idp1};
            Quadrant q2 (QuadBD, maxDepth, (*CurrentQuadIndex)++);
            q2.setIsTemplate(true);
            q2.setTemplateEdge(e1);

            Quadrants->push_back(q2);
            Quadrant *ptrQ2 = &Quadrants->back();
            TemplatesInsertados++;
            ti.quad2 = ptrQ2; 

            TemplesAdded.push_back(ti);
        }
    }


    bool Archipielago::SubSampling(const QuadEdge edge){
        unsigned int count = 0;
        Point3D p0 = (BackPoints)[edge[0]].getPoint();
        Point3D p1 = (BackPoints)[edge[1]].getPoint();

        for (const Quadrant &q : *Quadrants) {
            if (q.pointInside(*points, p0) && q.pointInside(*points, p1)) {
                return false;
            }
        }

        unsigned int samples = 10;
        for (unsigned int i = 1; i <= samples; ++i) {

            double t = double(i)/(samples+1);
            Point3D p = p0*(1.0-t) + p1*t;
            if (ply->pointIsInMesh(p)) {
                count++;       
            }
        }
        if (count>=1){
            return true;
        } 
        return false;
    }
    
    void Archipielago::getOutsideEdges(unsigned int maxDepth){

        for (const auto &entry : BackMapEdges) {
            const QuadEdge edge = entry.first;

            if (MapEdges.find(edge) != MapEdges.end())
                continue; 
            bool createTemplate = SubSampling(edge);

            if(createTemplate)createTemplates(edge, maxDepth);
        }
    }

bool Archipielago::pointInsideTemplate(const Quadrant &q, const Point3D &P) const {

        const auto &idx = q.getPointIndex();

        double minX = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double minY = std::numeric_limits<double>::max();
        double maxY = std::numeric_limits<double>::lowest();

        for (unsigned int i = 0; i < 4; ++i) {
            const Point3D &p = (*points)[idx[i]].getPoint();

            if (p.X() < minX) minX = p.X();
            if (p.X() > maxX) maxX = p.X();

            if (p.Y() < minY) minY = p.Y();
            if (p.Y() > maxY) maxY = p.Y();
        }

        if (P.X() > minX && P.X() < maxX &&
            P.Y() > minY && P.Y() < maxY) {
            return true;
        }

        return false;
    }

    QuadEdge Archipielago::SelectEdge(QuadEdge Tedge,QuadEdge Qedge1,QuadEdge Qedge2,unsigned int shared){

        auto getDirection = [&](const QuadEdge &e) -> Point3D {

            unsigned int other;

            if (e[0] == shared)
                other = e[1];
            else
                other = e[0];

            Point3D pShared = (*points)[shared].getPoint();
            Point3D pOther  = (*points)[other].getPoint();

            Point3D dir = pOther - pShared;
            dir.normalize();

            return dir;
        };

        Point3D dirT  = getDirection(Tedge);
        Point3D dirQ1 = getDirection(Qedge1);
        Point3D dirQ2 = getDirection(Qedge2);

        double dot1 = std::abs(dirT * dirQ1);
        double dot2 = std::abs(dirT * dirQ2);

        if (dot1 < dot2)
            return Qedge1;
        else
            return Qedge2;
    }

    unsigned int Archipielago::getIndexWithSharedVertex(Quadrant q,unsigned int shared,QuadEdge edge){

        vector<unsigned int> id = q.getPointIndex();
        unsigned int len = id.size();
        for(unsigned int i = 0; i < len; i++){
            if (id[i]== shared){
                if((id[i] == edge[0] && id[(i+1)%len] == edge[1]) || (id[i] == edge[1] && id[(i+1)%len] == edge[0])){
                    return id[(i-1)%len];
                }
                return id[(i+1)%len];
            }
        }
    }

    void Archipielago::fixTemplates(unsigned int maxDepth){
        ToFixWithQuads.clear();
        ToFixWithTems.clear();
   
        getTemplatesToFix();
        QuadAliasing qa;
        QuadEdge Qedge1,Qedge2,Quadrant_edge,Tedge;
        unsigned int sharedVertex;
        for(const FixWithQuad &ToRepair: ToFixWithQuads){
            sharedVertex = ToRepair.sharedVertex;
            Tedge = ToRepair.qTemplate.edge;
            const Quadrant &quad = *ToRepair.qQuadrant;
            qa.getAdjacentEdges(quad,sharedVertex,Qedge1,Qedge2);
            Quadrant_edge = SelectEdge(Tedge,Qedge1,Qedge2,sharedVertex);
            createFixTempQuad(Quadrant_edge,maxDepth,quad);
        }
        QuadEdge Tedge1, Tedge2;
        QuadEdge Edge1a,Edge1b,Edge2a,Edge2b;
        unsigned int p1a,p1b,p2a,p2b;

        for(FixWithTemp &ToRepair: ToFixWithTems){
            // cout << "Iteracion 1" <<endl;
            Tedge1 = ToRepair.Template1.edge;
            Tedge2 = ToRepair.Template2.edge;
            sharedVertex = ToRepair.sharedVertex;
            // cout << "Iteracion 2" <<endl;

            if (ToRepair.Template1.quad1) {
                const Quadrant quad1a = *ToRepair.Template1.quad1;
                p1a = getIndexWithSharedVertex(quad1a,sharedVertex,Tedge1);
            }
            if (ToRepair.Template1.quad2) {
                const Quadrant quad1b = *ToRepair.Template1.quad2;
                p1b = getIndexWithSharedVertex(quad1b,sharedVertex,Tedge1);
            }
            if (ToRepair.Template2.quad1) {
                const Quadrant quad2a = *ToRepair.Template2.quad1;
                p2a = getIndexWithSharedVertex(quad2a,sharedVertex,Tedge2);
            }

            if (ToRepair.Template2.quad2) {
                const Quadrant quad2b = *ToRepair.Template2.quad2;
                p2b = getIndexWithSharedVertex(quad2b,sharedVertex,Tedge2);
            }
            createFixWithTemp(sharedVertex,p1a,p1b,p2a,p2b,maxDepth);

        }

        cout<< "Templates Archipielagos Insertados: " << TemplatesInsertados << endl;
        ToFixWithQuads.clear();
        ToFixWithTems.clear();
        TemplesAdded.clear();



    }
}