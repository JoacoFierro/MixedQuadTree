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
#include "QuadAliasing.h"

 namespace Clobscode {
  
    void QuadAliasing::setQuadrant( list<Quadrant> &quadrants){
        this->Quadrants = &quadrants;

    }

    void QuadAliasing::setPoints( vector<MeshPoint> &Meshpoints){
        this->points = &Meshpoints;
    }

    void QuadAliasing::setActualIndex(unsigned int &idx){
        this->CurrentQuadIndex = &idx;
    }
    
    void QuadAliasing::QuadrantMap(){
        mapQuadrant.clear();

        for (const Quadrant &q : *Quadrants) {

            const auto &p = q.getPointIndex();

            for (size_t i=0;i<4;i++) {
                QuadEdge e(p[i], p[(i+1)%4]);
                mapQuadrant[e].push_back(&q);
            }
        }
    }



    void QuadAliasing::getAdjacentEdges(const Quadrant &q, unsigned int vertex, QuadEdge &e1, QuadEdge &e2){

        const auto &pts = q.getPointIndex();

        for (unsigned int i = 0; i < pts.size(); i++) {

            if (pts[i] != vertex)
                continue;

            unsigned int prev = (i + 3) % 4;
            unsigned int next = (i + 1) % 4;

            e1 = QuadEdge(pts[prev], pts[i]);
            e2 = QuadEdge(pts[i], pts[next]);

            return;
        }
        throw std::runtime_error("Shared vertex not found.");
    }

    bool  QuadAliasing::edgeHasNeighbour(const QuadEdge &e){
        auto it = mapQuadrant.find(e);

        if (it == mapQuadrant.end())
            return false;

        return it->second.size() > 1;
    }

    unsigned int QuadAliasing::oppositeVertex(const QuadEdge &e, unsigned int shared){
        return e[0]==shared ? e[1] : e[0];
    }
    
    bool QuadAliasing::ComparePoints(const Point3D &p1,const Point3D &p2){
            if (p1.X() == p2.X() && p1.Y() == p2.Y() && p1.Z() == p2.Z() )
                return true;
        return false;

    }
    void QuadAliasing::createTemplate(const VertexAlias &pinch){
        std::cout << "--FLAG 3----" << endl;
        unsigned int A = oppositeVertex(pinch.q1_edges[0],pinch.sharedVertex);
        unsigned int B = oppositeVertex(pinch.q1_edges[1],pinch.sharedVertex);
        unsigned int C = oppositeVertex(pinch.q2_edges[0],pinch.sharedVertex);
        unsigned int D = oppositeVertex(pinch.q2_edges[1],pinch.sharedVertex);

        Point3D PO = (*points)[pinch.sharedVertex].getPoint();
        Point3D PA = (*points)[A].getPoint();
        Point3D PB = (*points)[B].getPoint();
        Point3D PC = (*points)[C].getPoint();
        Point3D PD = (*points)[D].getPoint();
        Point3D dirOA = PA - PO;
        Point3D dirOB = PB - PO;
        Point3D dirOC = PC - PO;
        Point3D dirOD = PD - PO;

        dirOA.normalize();
        dirOB.normalize();
        dirOC.normalize();
        dirOD.normalize();
        Point3D nOA(-dirOA.Y(), dirOA.X(),0);
        Point3D nOB(-dirOB.Y(), dirOB.X(),0);
        Point3D nOC(-dirOC.Y(), dirOC.X(),0);
        Point3D nOD(-dirOD.Y(), dirOD.X(),0);

        nOA.normalize();
        nOB.normalize();
        nOC.normalize();
        nOD.normalize();
        double L =
            (PO-PA).Norm()+
            (PO-PB).Norm()+
            (PO-PC).Norm()+
            (PO-PD).Norm();

        L *= 0.25;

        double width = 0.45*L;

        //-----------------------------------------

        Point3D center1(0.0,0.0,0.0);
        Point3D center2(0.0,0.0,0.0); 

        const vector<unsigned int> ids1 = pinch.q1->getPointIndex(); 
        const vector<unsigned int> ids2 = pinch.q2->getPointIndex(); 

        for(unsigned int id : ids1){
            center1 = center1 + (*points)[id].getPoint();
        }
        std::cout << "FLAG 3.7" << endl;
        std::cout << "ID Cuadrante 2: " << const_cast<Clobscode::Quadrant*>(pinch.q2)->getIndex() <<endl;
        for(unsigned int id : ids2){
            std::cout << "FLAG 3.7.2.1" << endl;
            std::cout << id <<endl;
            center2 = center2 + (*points)[id].getPoint();
        }
        std::cout << "FLAG 3.8" << endl;
        center1 = center1 * 0.25;
        center2 = center2 * 0.25;
        
        //-----------------------------------------

        Point3D midpointOA=(PO+PA)*0.5;
        Point3D midpointOB=(PO+PB)*0.5;
        Point3D midpointOC=(PO+PC)*0.5;
        Point3D midpointOD=(PO+PD)*0.5;

        Point3D toCenterOA=center1-midpointOA;
        Point3D toCenterOB=center1-midpointOB;
        Point3D toCenterOC=center2-midpointOC;
        Point3D toCenterOD=center2-midpointOD;

        if(toCenterOA*nOA>0)nOA=nOA*(-1.0);
        if(toCenterOB*nOB>0)nOB=nOB*(-1.0);
        if(toCenterOC*nOC>0)nOC=nOC*(-1.0);
        if(toCenterOD*nOD>0)nOD=nOD*(-1.0);

        Point3D bevelAIn = dirOA + nOA;
        Point3D bevelBIn = dirOB + nOB;
        Point3D bevelCIn = dirOC + nOC;
        Point3D bevelDIn = dirOD + nOD;

        Point3D bevelAOut = -dirOA + nOA;
        Point3D bevelBOut = -dirOB + nOB;
        Point3D bevelCOut = -dirOC + nOC;
        Point3D bevelDOut = -dirOD + nOD;

        bevelAIn.normalize();
        bevelBIn.normalize();
        bevelCIn.normalize();
        bevelDIn.normalize();

        bevelAOut.normalize();
        bevelBOut.normalize();
        bevelCOut.normalize();
        bevelDOut.normalize();

        
        Point3D Aout = PA + bevelAOut*width;
        Point3D Bout = PB + bevelBOut*width;
        Point3D Cout = PC + bevelCOut*width;
        Point3D Dout = PD + bevelDOut*width;

        Point3D Ai = PO + bevelAIn*width;
        Point3D Bi = PO + bevelBIn*width;
        Point3D Ci = PO + bevelCIn*width;
        Point3D Di = PO + bevelDIn*width;

        //--------------------------Implementar descarte de templates si es que existen


        //------------------------ Im

        unsigned int idAout=points->size();
        points->push_back(MeshPoint(Aout));

        unsigned int idBout=points->size();
        points->push_back(MeshPoint(Bout));

        unsigned int idCout=points->size();
        points->push_back(MeshPoint(Cout));

        unsigned int idDout=points->size();
        points->push_back(MeshPoint(Dout));


        Point3D CommonPointA = ComparePoints(Ai,Ci) ? Ci : Di;
        Point3D CommonPointB = ComparePoints(Ai,Ci) ? Di : Ci;
        
        unsigned int idCommonA=points->size();
        points->push_back(MeshPoint(CommonPointA));

        unsigned int idCommonB=points->size();
        points->push_back(MeshPoint(CommonPointB));


        vector<unsigned int> QuadA = {pinch.sharedVertex,A,idAout,idCommonA};
        vector<unsigned int> QuadB = {pinch.sharedVertex,B,idBout,idCommonB};
        vector<unsigned int> QuadC = {pinch.sharedVertex,C,idCout, ComparePoints(Ai,Ci) ? idCommonA: idCommonB };
        vector<unsigned int> QuadD = {pinch.sharedVertex,D,idDout, ComparePoints(Ai,Ci) ? idCommonB: idCommonA };

        unsigned short ref_level1 = pinch.q1->getRefinementLevel();
        unsigned short ref_level2 = pinch.q2->getRefinementLevel();

        Quadrant qa (QuadA, ref_level1, (*CurrentQuadIndex)++);
        Quadrants->push_back(qa);

        Quadrant qb (QuadB, ref_level1, (*CurrentQuadIndex)++);
        Quadrants->push_back(qb);

        Quadrant qc (QuadC, ref_level2, (*CurrentQuadIndex)++);
        Quadrants->push_back(qc);

        Quadrant qd (QuadD, ref_level2, (*CurrentQuadIndex)++);
        Quadrants->push_back(qd);
        std::cout << "=================================" << endl;
    }   

    void QuadAliasing::CreateTemplates(){
        for(VertexAlias pinch: Pinches){
            createTemplate(pinch);
        }
    }
    
    void QuadAliasing::getPinches(){
        QuadrantMap();
        list<Quadrant>::iterator iter1;
        list<Quadrant>::iterator iter2;
        for (iter1 = Quadrants->begin(); iter1 != Quadrants->end(); ++iter1){
            
            auto iter2 = iter1;
            ++iter2;

            for (; iter2 != Quadrants->end(); ++iter2){
                int common = 0;
                unsigned int sharedVertex;

                for (unsigned int a : iter1->getPointIndex()){
                    for (unsigned int b : iter2->getPointIndex()){
                        if (a == b){
                            common++;
                            sharedVertex = a;
                        }
                    }
                }
                if(common == 1){
                    QuadEdge e1,e2;
                    getAdjacentEdges(*iter1, sharedVertex,e1,e2);

                    bool flag1 = edgeHasNeighbour(e1);
                    bool flag2 = edgeHasNeighbour(e2);
                    if (flag1 || flag2){
                        continue;
                    }
                    VertexAlias Pinch;
                    Pinch.sharedVertex = sharedVertex;
                    Pinch.q1 = &(*iter1);
                    Pinch.q2 = &(*iter2);
                    Pinch.q1_edges = {e1,e2};
                    getAdjacentEdges(*iter2, sharedVertex, Pinch.q2_edges[0], Pinch.q2_edges[1]);
                    Pinches.push_back(Pinch);
                }
                else{
                    continue;
                }
            }
        }
    }

    void QuadAliasing::printPinches(){
        std::cout << "\n========== PINCHES DETECTADOS ==========\n";

        int idx = 0;
        for (const VertexAlias &pinche : Pinches){
        
            std::cout << "\nPinch " << idx++ << "\n";
            std::cout << "----------------------------------------\n";

            std::cout << "Shared vertex: "
                    << pinche.sharedVertex << "\n";

            std::cout << "Quadrant 1: ";
            for (unsigned int p : pinche.q1->getPointIndex())
                std::cout << p << " ";
            std::cout << "\n";

            std::cout << "Quadrant 2: ";
            for (unsigned int p : pinche.q2->getPointIndex())
                std::cout << p << " ";
            std::cout << "\n";

            std::cout << "Adjacent edges of Q1:\n";
            std::cout << "    " << pinche.q1_edges[0] << "\n";
            std::cout << "    " << pinche.q1_edges[1] << "\n";

            std::cout << "Adjacent edges of Q2:\n";
            std::cout << "    " << pinche.q2_edges[0] << "\n";
            std::cout << "    " << pinche.q2_edges[1] << "\n";
        }
        std::cout << "\nTotal pinches: " << Pinches.size() << "\n";
        std::cout << "========================================\n";
    }
    void QuadAliasing::savePinches(const std::string &filename){

        std::ofstream file(filename);

        if (!file.is_open()) {
            std::cerr << "Error: no se pudo crear el archivo "
                    << filename << std::endl;
            return;
        }

        file << "========== PINCHES DETECTADOS ==========\n";

        int idx = 0;

        for (const VertexAlias &pinche : Pinches) {

            file << "\nPinch " << idx++ << "\n";
            file << "----------------------------------------\n";

            file << "Shared vertex: "
                << pinche.sharedVertex << "\n";

            file << "Quadrant 1: ";
            for (unsigned int p : pinche.q1->getPointIndex())
                file << p << " ";
            file << "\n";

            file << "Quadrant 2: ";
            for (unsigned int p : pinche.q2->getPointIndex())
                file << p << " ";
            file << "\n";

            file << "Index Quadrant 2: ";
            file << const_cast<Clobscode::Quadrant*>(pinche.q2)->getIndex() << "\n";

            file << "Adjacent edges of Q1:\n";
            file << "    " << pinche.q1_edges[0] << "\n";
            file << "    " << pinche.q1_edges[1] << "\n";

            file << "Adjacent edges of Q2:\n";
            file << "    " << pinche.q2_edges[0] << "\n";
            file << "    " << pinche.q2_edges[1] << "\n";
        }

        file << "\nTotal pinches: " << Pinches.size() << "\n";
        file << "========================================\n";

        file.close();

        std::cout << "Pinches guardados en: "
                << filename << std::endl;
    }


    void QuadAliasing::printPoints(){
        std::cout << "========== Points ======== \n";
        int idx = 0;
        for (const MeshPoint p: *points){
            std::cout << "-----Point "<< idx++ <<" : " 
                      << p.getPoint() << "\n";

            std::cout << "Elements:";
            for(unsigned int e: p.getElements()){
                std::cout << e << " ";
            }
            std::cout << "\n";
            std::cout << "MaxDistance:"
                      << p.getMaxDistance() << "\n";      
        }
    }

void QuadAliasing::printQuadrants()
{
    std::cout << "========== Quadrants ==========\n";

    int idx = 0;

    for (Quadrant &q : *Quadrants)
    {
        std::cout << "\nQuadrant " << idx++ << '\n';
        std::cout << "--------------------------------\n";

        // Corner points
        std::cout << "Corner points: ";
        for (unsigned int p : q.getPointIndex())
            std::cout << p << " ";
        std::cout << "\n";

        std::cout << "Quadrant index: "
                  << q.getIndex() << "\n";

        // Refinement level
        std::cout << "Refinement level: "
                  << q.getRefinementLevel() << "\n";

        // Inside / Surface
        std::cout << "Inside: "
                  << std::boolalpha << q.isInside() << "\n";

        std::cout << "Intersects surface: "
                  << std::boolalpha << q.intersectsSurface() << "\n";

        std::cout << "Surface: "
                  << std::boolalpha << q.isSurface() << "\n";

        std::cout << "Debugging: "
                  << std::boolalpha << q.isDebugging() << "\n";

        // Max distance
        std::cout << "Max distance: "
                  << q.getMaxDistance() << "\n";

        // Intersected edges
        std::cout << "Intersected edges (" 
                  << q.getIntersectedEdges().size() << "): ";

        for (unsigned int e : q.getIntersectedEdges())
            std::cout << e << " ";

        std::cout << "\n";

        // Intersected features
        std::cout << "Intersected features ("
                  << q.getIntersectedFeatures().size() << "): ";

        for (unsigned int f : q.getIntersectedFeatures())
            std::cout << f << " ";

        std::cout << "\n";

        // All subpoints
        std::cout << "Sub points: ";

        for (unsigned int p : q.getSubPointIndex())
            std::cout << p << " ";

        std::cout << "\n";

        // Edge subpoints
        std::cout << "Sub points on edges: ";

        for (unsigned int p : q.getSubPointIndexOnlyOnEdges())
            std::cout << p << " ";

        std::cout << "\n";

        // Sub elements
        const auto &subs = q.getSubElements();

        std::cout << "Sub elements (" << subs.size() << ")\n";

        for (size_t i = 0; i < subs.size(); ++i)
        {
            std::cout << "  [" << i << "] ";

            for (unsigned int p : subs[i])
                std::cout << p << " ";

            std::cout << "\n";
        }
    }

    std::cout << "Cantidad de cuadrantes: " << Quadrants->size() <<"\n";

    std::cout << "===============================\n";
}
}


 