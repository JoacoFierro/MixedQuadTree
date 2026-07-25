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

    void QuadAliasing::setInput(const Polyline &ply){
        this->ply = &ply;
    }

    void QuadAliasing::setActualIndex(unsigned int idx){
        this->CurrentQuadIndex = idx;
    }
    
    void QuadAliasing::QuadrantVertexMap(){
        mapQuadrant.clear();
        mapVertexQuadrants.clear();

        for (const Quadrant &q : *Quadrants) {

            const auto &p = q.getPointIndex();

            for (unsigned int i=0;i< p.size();i++) {
                QuadEdge e(p[i], p[(i+1)%4]);
                mapQuadrant[e].push_back(&q);
                mapVertexQuadrants[p[i]].push_back(&q);
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
        constexpr double eps = 1e-6;

        return std::abs(p1.X() - p2.X()) < eps &&
            std::abs(p1.Y() - p2.Y()) < eps &&
            std::abs(p1.Z() - p2.Z()) < eps;

    }


    void QuadAliasing::createTemplate(const VertexAlias &pinch){
        
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
        for(unsigned int id : ids2){
            center2 = center2 + (*points)[id].getPoint();
        }

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

        Point3D CommonPointA = ComparePoints(Ai,Ci) ? Ci : Di;
        Point3D CommonPointB = ComparePoints(Ai,Ci) ? Di : Ci;
        
        bool InsertA = true, InsertB = true,InsertC = true , InsertD = true;
        bool InsertAinVector = true, InsertBinVector = true,InsertCinVector = true, InsertDinVector = true, 
                                                            InsertCommonAV = true,InsertCommonBV = true; 
        bool InsertCommonA = true , InsertCommonB = true;

        unsigned int idAout,idBout,idCout,idDout;
        unsigned int idCommonA = 0;
        unsigned int idCommonB = 0;

        for (unsigned int i = 0; i < points->size(); ++i) {
            const Point3D &p = (*points)[i].getPoint();
            if (ComparePoints(Aout,p)){
                idAout = i;
                InsertAinVector = false;
            }
            if (ComparePoints(Bout,p)){
                idBout = i;
                InsertBinVector = false;
            }
            if (ComparePoints(Cout,p)){
                idCout = i;
                InsertCinVector = false;
            }
            if (ComparePoints(Dout,p)){
                idDout = i;
                InsertDinVector = false;
            }
            if (ComparePoints(CommonPointA,p)){
                idCommonA = i;
                InsertCommonAV = false;
            }
            if (ComparePoints(CommonPointB,p)){
                idCommonB = i;
                InsertCommonBV = false;
            }
        }
        vector<unsigned int> indices,indA,indB,indC,indD;

        for (const Quadrant &q : *Quadrants) {
            if (q.getIsTemplate()){
   
                vector<unsigned int> indices =q.getPointIndex();
                indA = {pinch.sharedVertex,A,idAout,idCommonA};
                indB = {pinch.sharedVertex,B,idBout,idCommonB};
                indC = {pinch.sharedVertex,C,idCout, ComparePoints(Ai,Ci) ? idCommonA: idCommonB};
                indD = {pinch.sharedVertex,D,idDout, ComparePoints(Ai,Ci) ? idCommonB: idCommonA};

                std::sort(indices.begin(),indices.end());
                std::sort(indA.begin(),indA.end());
                std::sort(indB.begin(),indB.end());
                std::sort(indC.begin(),indC.end());
                std::sort(indD.begin(),indD.end());

                if(indices == indA)InsertA = false;
                if(indices == indB)InsertB = false;
                if(indices == indC)InsertC = false;
                if(indices == indD)InsertD = false;

                if (!InsertA  && !(ComparePoints(Ai,Ci) ? InsertC: InsertD))InsertCommonA = false;
                if (!InsertB  && !(ComparePoints(Ai,Ci) ? InsertD: InsertC))InsertCommonB = false;

                continue;
            }
            if (q.pointInside(*points, Aout)) {
                InsertA = false;
            }
            if (q.pointInside(*points, Bout)) {
                InsertB = false;
            }
            if (q.pointInside(*points, Cout)) {
                InsertC = false;
            }
            if (q.pointInside(*points, Dout)) {
                InsertD = false;
            }
            if (q.pointInside(*points, CommonPointA)) {
                InsertCommonA = false;
            }
            if (q.pointInside(*points, CommonPointB)) {
                InsertCommonB = false;
            }
        }

        if(InsertCommonAV && InsertCommonA){
            idCommonA=points->size();
            points->push_back(MeshPoint(CommonPointA));
        }
        if(InsertCommonBV && InsertCommonB){
           idCommonB=points->size();
           points->push_back(MeshPoint(CommonPointB));
        }
        if(InsertAinVector && InsertA){
            idAout=points->size();
            points->push_back(MeshPoint(Aout));
        }
        if(InsertBinVector && InsertB){
            idBout=points->size();
            points->push_back(MeshPoint(Bout));
        }
        if(InsertCinVector && InsertC){
            idCout=points->size();
            points->push_back(MeshPoint(Cout));
        }
        if(InsertDinVector && InsertD){
            idDout=points->size();
            points->push_back(MeshPoint(Dout));
        }

        unsigned short ref_level1 = 0;
        unsigned short ref_level2 = 0;

        if(InsertA || InsertB){
            ref_level1 = pinch.q1->getRefinementLevel();
        }
        if(InsertC || InsertD){
            ref_level2 = pinch.q2->getRefinementLevel();
        }

        if(InsertA && InsertCommonA){
            vector<unsigned int> QuadA = {pinch.sharedVertex,A,idAout,idCommonA};
            Quadrant qa (QuadA, ref_level1, CurrentQuadIndex++);
            qa.setIsTemplate(true);
            Quadrants->push_back(qa);
        }
        if(InsertB && InsertCommonB){
            vector<unsigned int> QuadB = {pinch.sharedVertex,B,idBout,idCommonB};
            Quadrant qb (QuadB, ref_level1, CurrentQuadIndex++);
            qb.setIsTemplate(true);
            Quadrants->push_back(qb);
            
        }
        if(InsertC && (ComparePoints(Ai,Ci) ? InsertCommonA: InsertCommonB)){
            vector<unsigned int> QuadC = {pinch.sharedVertex,C,idCout, ComparePoints(Ai,Ci) ? idCommonA: idCommonB };
            Quadrant qc (QuadC, ref_level2, CurrentQuadIndex++);
            qc.setIsTemplate(true);
            Quadrants->push_back(qc);
        }
        if(InsertD && (ComparePoints(Ai,Ci) ? InsertCommonB: InsertCommonA)){
            vector<unsigned int> QuadD = {pinch.sharedVertex,D,idDout, ComparePoints(Ai,Ci) ? idCommonB: idCommonA };
            Quadrant qd (QuadD, ref_level2, CurrentQuadIndex++);
            qd.setIsTemplate(true);
            Quadrants->push_back(qd);
        }

}   

    void QuadAliasing::CreateTemplates(){
        for(VertexAlias pinch: Pinches){
            createTemplate(pinch);
        }
    }
    
    void QuadAliasing::getPinches(){
        QuadrantVertexMap();
        cout<<"Pinches Encontrados: "<< mapVertexQuadrants.size() <<endl;
        for (auto &entry : mapVertexQuadrants){
        
            unsigned int sharedVertex = entry.first;
            auto &incident = entry.second;

            if (incident.size() != 2)
                continue;
        
            const Quadrant *q1 = incident[0];
            const Quadrant *q2 = incident[1];
            QuadEdge q1e1,q1e2;
            QuadEdge q2e1,q2e2;

            getAdjacentEdges(*q1,sharedVertex,q1e1,q1e2);
            getAdjacentEdges(*q2,sharedVertex,q2e1,q2e2);

            if(edgeHasNeighbour(q1e1)) continue;
            if(edgeHasNeighbour(q1e2)) continue;
            if(edgeHasNeighbour(q2e1)) continue;
            if(edgeHasNeighbour(q2e2)) continue;

            unsigned int A = oppositeVertex(q1e1,sharedVertex);
            unsigned int B = oppositeVertex(q1e2,sharedVertex);
            unsigned int C = oppositeVertex(q2e1,sharedVertex);
            unsigned int D = oppositeVertex(q2e2,sharedVertex);

            if(!validAngularConfiguration(sharedVertex,A,B,C,D))
                continue;

            VertexAlias pinch;
            pinch.sharedVertex=sharedVertex;
            pinch.q1=q1;
            pinch.q2=q2;
            pinch.q1_edges={q1e1,q1e2};
            pinch.q2_edges={q2e1,q2e2};
            Pinches.push_back(pinch);
        }
    }

    bool QuadAliasing::validAngularConfiguration(
            unsigned int shared,
            unsigned int A,
            unsigned int B,
            unsigned int C,
            unsigned int D)
    {
        const Point3D &O = (*points)[shared].getPoint();

        std::vector<double> ang;
        ang.reserve(4);

        unsigned int ids[4] = {A,B,C,D};

        for (unsigned int id : ids) {

            Point3D v = (*points)[id].getPoint() - O;

            double a = atan2(v.Y(), v.X());

            if (a < 0.0)
                a += 2.0*M_PI;

            ang.push_back(a);
        }

        std::sort(ang.begin(), ang.end());

        double maxGap = 0.0;

        for (int i=0;i<4;i++) {

            double next = ang[(i+1)%4];

            if (i==3)
                next += 2.0*M_PI;

            double gap = next - ang[i];

            if (gap > maxGap)
                maxGap = gap;
        }

        // Si existe un hueco muy grande,
        // los cuatro vecinos NO rodean al vértice.

        return maxGap < M_PI;      // 180°
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

        file << "========== PINCHES CON DISTINTO NIVEL DE REFINAMIENTO ==========\n";

        int idx = 0;

        for (const VertexAlias &pinche : Pinches) {

            unsigned short level1 = pinche.q1->getRefinementLevel();
            unsigned short level2 = pinche.q2->getRefinementLevel();

            // Solo guardar los casos donde los niveles son distintos
            if (level1 == level2)
                continue;

            file << "\nPinch " << idx++ << "\n";
            file << "--------------------------------------------------\n";

            file << "Shared vertex: " << pinche.sharedVertex << "\n\n";

            file << "Quadrant 1\n";
            file << "  Index: " << const_cast<Clobscode::Quadrant*>(pinche.q1)->getIndex() << "\n";
            file << "  Refinement level: " << level1 << "\n";
            file << "  Vertices: ";
            for (unsigned int p : pinche.q1->getPointIndex())
                file << p << " ";
            file << "\n";

            file << "  Adjacent edges:\n";
            file << "    " << pinche.q1_edges[0] << "\n";
            file << "    " << pinche.q1_edges[1] << "\n\n";

            file << "Quadrant 2\n";
            file << "  Index: " << const_cast<Clobscode::Quadrant*>(pinche.q2)->getIndex() << "\n";
            file << "  Refinement level: " << level2 << "\n";
            file << "  Vertices: ";
            for (unsigned int p : pinche.q2->getPointIndex())
                file << p << " ";
            file << "\n";

            file << "  Adjacent edges:\n";
            file << "    " << pinche.q2_edges[0] << "\n";
            file << "    " << pinche.q2_edges[1] << "\n";

            file << "--------------------------------------------------\n";
            unsigned int A = oppositeVertex(pinche.q1_edges[0], pinche.sharedVertex);
            unsigned int B = oppositeVertex(pinche.q1_edges[1], pinche.sharedVertex);
            unsigned int C = oppositeVertex(pinche.q2_edges[0], pinche.sharedVertex);
            unsigned int D = oppositeVertex(pinche.q2_edges[1], pinche.sharedVertex);

            file << "Opposite vertices:\n";
            file << "  A = " << A << "\n";
            file << "  B = " << B << "\n";
            file << "  C = " << C << "\n";
            file << "  D = " << D << "\n";
        }

        file << "\nTotal pinches escritos: " << idx << "\n";
        file << "==================================================\n";

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


 