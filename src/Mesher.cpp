/*
 <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>
 
 Copyright (C) <2013,2020>  <Claudio Lobos> All rights reserved.
 
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
 * @file Mesher.cpp
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief
 **/

#include "Mesher.h"
#include <math.h>
#include <iomanip>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <functional>

namespace Clobscode
{
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::Mesher() {}
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    Mesher::~Mesher() {}
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 4^rl times (where rl stands
    //for Refinement Level).
    std::shared_ptr<FEMesh> Mesher::refineMesh(Polyline &input, const unsigned short &rl,
                                               const string &name, list<unsigned int> &roctli,
                                               list<RefinementRegion *> &all_reg,
                                               GeometricTransform &gt,
                                               const unsigned short &minrl,
                                               const unsigned short &omaxrl,
                                               const bool &debugging, unsigned int sampleSize,
                                               bool decoration,
                                               bool useTusqh,
                                               unsigned int tusqhSampleSize,
                                               bool refineOnEdgeIntersect,
                                               unsigned int tusqhExtraResolveDepth,bool Aliasing,
                                               bool useSubgrid,
                                               unsigned int subgridSampleSize,
                                               double subgridJoinThreshold,
                                               unsigned int subgridMinComponentCells){

        //Note: rotation are not enabled when refining an already produced mesh.
        bool rotated = !gt.Default();
        if(rotated) {
            /*cout << "rotating input surface mesh\n";
             cout << gt.getCentroid() << "\n";
             cout << gt.getXAxis() << " " << gt.getYAxis();
             cout << " " << gt.getZAxis() << "\n";*/
            gt.rotatePolyline(input);
        }

#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> grid_octree=make_shared<FEMesh>();
            saveOutputMesh(grid_octree,points,Quadrants,debugging);
            string tmp_name = name + "_grid";
            Services::WriteVTK(tmp_name,grid_octree);
        }
#endif

        //split Quadrants until the refinement level (rl) is achieved.
        //The output will be a one-irregular mesh.
        if (useTusqh) {
            // Heuristic for the TUSQH coarse-grid bug: see Mesher::
            // generateMesh. Here the quadtree was loaded from a
            // previous mesh (option -c), so the initial grid may be
            // very coarse if that mesh was produced with `generate
            // Mesh` and -T but without enough subdivisions (e.g. the
            // Chesapeake Bay case). Apply the same heuristic.
            preRefineForTusqh(input, rl, all_reg, name, debugging, Aliasing);
            windingSubdivide(input, rl, tusqhSampleSize,
                             refineOnEdgeIntersect, name, debugging,
                             tusqhExtraResolveDepth);
        } else {
            splitQuadrants(rl,input,roctli,all_reg,name,minrl,omaxrl,debugging);
        }

        // compute volume fractions using winding numbers with s x s samples
        mSampleSize = sampleSize;
        computeVolumeFractions(input, mSampleSize);

        // TUSQH §3.3 sub-cell volume fractions + §3.4 archipelago resolution.
        // Both must run BEFORE removing surface quads because archipelago
        // detection needs the full quad connectivity.
        if (useSubgrid) {
            computeSubcellVolumeFractions(input, subgridSampleSize,
                                           subgridJoinThreshold, name);
            resolveArchipelagos(input, subgridSampleSize,
                                subgridJoinThreshold,
                                subgridMinComponentCells,
                                name,Aliasing);
        }

        //Save the Octant mesh for further refinement.
        Services::WriteQuadtreeMesh(name,points,Quadrants,MapEdges,gt);
        //Some Quads will be then removed due to proximity with the surface.
        //However we must preserve them if the oct mesh to avoid congruency
        //problems. For this reason we will keep track of removed quads
        //so we can easily link elements to quad index when reading an oct
        //mesh file.
        map<unsigned int, bool> removedquads;
        list<unsigned int> quadmeshidx;
        for (auto q: Quadrants) {
            removedquads[q.getIndex()] = true;
            quadmeshidx.push_back(q.getIndex());
        }
        
        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);
        computeNodeMaxDist();
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree=make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants,debugging);
            string tmp_name = name + "_quads";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        projectCloseToBoundaryNodes(input);

#if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> closeto_octree=make_shared<FEMesh>();
            saveOutputMesh(closeto_octree,points,Quadrants);
            string tmp_name = name + "_closeto";
            Services::WriteVTK(tmp_name,closeto_octree);
        }
#endif
        
        removeOnSurfaceSafe(input);
        
        //update element and node info.
        linkElementsToNodes();
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree=make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_remSur";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        //shrink outside nodes to the input domain boundary
        //[DISABLED] shrinkToBoundary(input);
        
#if (VTKOUT==true)        //CL Debbuging
        //[DISABLED]
        //    std::shared_ptr<FEMesh> shrink_octree=make_shared<FEMesh>();
        //    saveOutputMesh(shrink_octree,points,Quadrants);
        //    string tmp_name = name + "_shrink";
        //    Services::WriteVTK(tmp_name,shrink_octree);
#endif

        //apply the surface Patterns
        //[DISABLED] applySurfacePatterns(input);
        //removeOnSurface(input);
        
        if (rotated) {
            // rotate the mesh
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
            // rotate back the polyline as well
            //FJA no need for Refinement ??
            //for (unsigned int i=0; i<input.getPoints().size(); i++) {
            //    gt.applyInverse(input.getPoints()[i]);
            //}
        }
        
        //the almighty output mesh
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh>();
        
        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        //Write element-quad info the file
        Services::addOctElemntInfo(name,Quadrants,removedquads,quadmeshidx);
        
        return mesh;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //Create a grid mesh regarding the Bounding Box of input mesh.
    //This will produce several cubes as roots of an octree structure.
    //Then split each initial element 4^rl times (where rl stands
    //for Refinement Level).
    std::shared_ptr<FEMesh> Mesher::generateMesh(Polyline &input, const unsigned short &rl,
                                                 const string &name,
                                                 list<RefinementRegion *> &all_reg,
                                                 const bool &debugging, unsigned int sampleSize,
                                                 bool decoration,
                                                 bool useTusqh,
                                                 unsigned int tusqhSampleSize,
                                                 bool refineOnEdgeIntersect,
                                                 unsigned int tusqhExtraResolveDepth, bool Aliasing,
                                                 bool useSubgrid,
                                                 unsigned int subgridSampleSize,
                                                 double subgridJoinThreshold,
                                                 unsigned int subgridMinComponentCells){

        //ATTENTION: geometric transform causes invalid input rotation when the
        //input is a cube.
        GeometricTransform gt;

        //rotate: This method is written below and its mostly commented because
        //it causes conflicts when the input is a cube. Must be checked.
        bool rotated = rotateGridMesh(input, all_reg, gt);

        //generate root Quadrants
        generateGridMesh(input);

#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> grid_octree = make_shared<FEMesh>();
            saveOutputMesh(grid_octree,points,Quadrants,debugging);
            string tmp_name = name + "_grid";
            Services::WriteVTK(tmp_name,grid_octree);
        }
#endif

        //split Quadrants until the refinement level (rl) is achieved.
        //the last 0 correspond to min RL in the mesh. As this mesh
        //is starting from scratch, this value is set to 0. When refining
        //an existing mesh, this value may change.
        //The output will be a one-irregular mesh.
        if (useTusqh) {
            // Heurística para el error de cuadrícula gruesa de TUSQH: si el quadtree
            // comienza con muy pocas celdas raíz (1 o 2 es lo típico para polilíneas
            // grandes como la Bahía de Chesapeake en el nivel 0) y el parámetro -N es pequeño
            // (por defecto 2 -> 4 muestras por celda), TUSQH submuestrea las
            // celdas enormes y las clasifica erróneamente como TodoAdentro/TodoAfuera (AllInside/AllOutside),
            // lo que produce casi nulas subdivisiones. Se pre-refina uniformemente hasta
            // un nivel base usando el método clásico `generateQuadtreeMesh` para que
            // TUSQH tenga suficientes celdas pequeñas y sus cálculos sean significativos.
            preRefineForTusqh(input, rl, all_reg, name, debugging, Aliasing);
            windingSubdivide(input, rl, tusqhSampleSize,
                             refineOnEdgeIntersect, name, debugging,
                             tusqhExtraResolveDepth,Aliasing);
        } else {
            generateQuadtreeMesh(rl,input,all_reg,name,0,debugging,Quadrants.size(),Aliasing);
        }

#if (VTKOUT==true)
        {
            // Step 2 — Post-TUSQH quadtree state (BEFORE volume fractions).
            // This snapshot shows how TUSQH subdivided the polyline, before
            // any sub-cell VF / archipelago work. It is the baseline that
            // the pre/post-archipelago VTKs (steps 4.5 and 5) are diffed
            // against.
            std::shared_ptr<FEMesh> octree_mesh = make_shared<FEMesh>();
            if (!Quadrants.empty()) {
                saveOutputMesh(octree_mesh, points, Quadrants);
            }
            string tmp_name = name + "_octree";
            Services::WriteVTK(tmp_name, octree_mesh);
        }
#endif

        // compute volume fractions using winding numbers with s x s samples
        mSampleSize = sampleSize;
        computeVolumeFractions(input, mSampleSize);

        // TUSQH §3.3 sub-cell volume fractions + §3.4 archipelago resolution.
        // Both must run BEFORE removing surface quads because archipelago
        // detection needs the full quad connectivity.
        if (useSubgrid) {
            computeSubcellVolumeFractions(input, subgridSampleSize,
                                           subgridJoinThreshold, name);

#if (VTKOUT==true)
            {
                // Step 4.5 — Pre-archipelago state. After sub-cell VF has
                // been computed but BEFORE bridges are added. The edge
                // sub-cell VF array (mEdgeSubcellVF) is what
                // SubgridSampler::buildQuadPerpThickness produced; if
                // that function ever reads `Quadrants[info[k]]` with a
                // q_id that is not equal to the vector index, the perp
                // thickness (and therefore the edge VF) will be wrong.
                // Diffing this file against _postarchipelago shows the
                // exact effect of the bridge loop.
                std::shared_ptr<FEMesh> pre_mesh = make_shared<FEMesh>();
                if (!Quadrants.empty()) {
                    saveOutputMesh(pre_mesh, points, Quadrants);
                }
                string tmp_name = name + "_prearchipelago";
                Services::WriteVTK(tmp_name, pre_mesh);
            }
#endif

            resolveArchipelagos(input, subgridSampleSize,
                                subgridJoinThreshold,
                                subgridMinComponentCells,
                                name,Aliasing);
        }

        Services::WriteQuadtreeMesh(name,points,Quadrants,MapEdges,gt);
        //Some Quads will be then removed due to proximity with the surface.
        //However we must preserve them if the oct mesh to avoid congruency
        //problems. For this reason we will keep track of removed quads
        //so we can easily link elements to quad index when reading an oct
        //mesh file.
        map<unsigned int, bool> removedquads;
        list<unsigned int> quadmeshidx;
        for (auto q: Quadrants) {
            removedquads[q.getIndex()] = true;
            quadmeshidx.push_back(q.getIndex());
        }
        
        
        //link element and node info for code optimization, also
        //detect Quadrants with features.
        detectFeatureQuadrants(input);
        linkElementsToNodes();
        detectInsideNodes(input);
        computeNodeMaxDist();
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree = make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_quads";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        projectCloseToBoundaryNodes(input);
        
#if (VTKOUT==true)         //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> closeto_octree = make_shared<FEMesh>();
            saveOutputMesh(closeto_octree,points,Quadrants);
            string tmp_name = name + "_closeto";
            Services::WriteVTK(tmp_name,closeto_octree);
        }
#endif

        removeOnSurfaceSafe(input);
        
        //update element and node info.
        linkElementsToNodes();
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> pure_octree = make_shared<FEMesh>();
            saveOutputMesh(pure_octree,points,Quadrants);
            string tmp_name = name + "_remSur";
            Services::WriteVTK(tmp_name,pure_octree);
        }
#endif

        //shrink outside nodes to the input domain boundary
        //[DISABLED] shrinkToBoundary(input);
        
#if (VTKOUT==true)        //CL Debbuging
        //[DISABLED]
        //    std::shared_ptr<FEMesh> shrink_octree = make_shared<FEMesh>();
        //    saveOutputMesh(shrink_octree,points,Quadrants);
        //    string tmp_name = name + "_shrink";
        //    Services::WriteVTK(tmp_name,shrink_octree);
#endif

        //apply the surface Patterns
        //[DISABLED] applySurfacePatterns(input);
        
        if (rotated) {
            // rotate the mesh
            for (unsigned int i=0; i<points.size(); i++) {
                gt.applyInverse(points[i].getPoint());
            }
            // rotate back the polyline as well
            for (unsigned int i=0; i<input.getPoints().size(); i++) {
                gt.applyInverse(input.getPoints()[i]);
            }
        }
        
        //the almighty output mesh
        std::shared_ptr<FEMesh> mesh = std::make_shared<FEMesh>();
        
        //save the data of the mesh in its final state
        saveOutputMesh(mesh,decoration);
        
        //Write element-quad info the file
        Services::addOctElemntInfo(name,Quadrants,removedquads,quadmeshidx);
        
        return mesh;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool Mesher::rotateGridMesh(Polyline &input, list<RefinementRegion *> &all_reg,
                                GeometricTransform &gt){
        
        list<RefinementRegion *>::const_iterator it, rrot;
        bool inputHasbeenRotated = false;
        
        for (it = all_reg.begin(); it!=all_reg.end(); it++) {
            //in case of input roi
            if((*it)->needsInputRotation()){
                if(!inputHasbeenRotated){
                    gt = (*it)->rotateWithinYou(input);
                    inputHasbeenRotated = true;
                    rrot=it;
                    break;
                }
            }
        }
        
        if (inputHasbeenRotated) {
            for (it = all_reg.begin(); it!=all_reg.end(); it++) {
                if (it!=rrot) {
                    if ((*it)->needsLocalRotation()) {
                        (*it)->rotate(gt);
                    }
                }
            }
        }
        
        return inputHasbeenRotated;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::generateGridMesh(Polyline &input){
        //vectors with each coordinate per axis
        vector<double> all_x, all_y;
        vector<vector<unsigned int> > elements;
        
        auto start_time = chrono::high_resolution_clock::now();
        
        GridMesher gm;
        gm.generatePoints(input.getBounds(),all_x,all_y);
        gm.generateMesh(all_x,all_y,points,elements);
        
        Quadrants.reserve(elements.size());
        
        //create the root Quadrants
        for (unsigned int i=0; i<elements.size(); i++) {
            Quadrant o (elements[i], 0, i);
            //Only when the Quadrant intersects the input
            //add it to the list of current Quadrants. As
            //This is the first time Quadrants are checked
            //for intersections they must be made w.r.t.
            //all input edges.
            IntersectionsVisitor iv(false);
            //if (o.checkIntersections(input,points)) {
            iv.setPolyline(input);
            iv.setPoints(points);
            if (o.accept(&iv)) {
                EdgeVisitor::insertEdges(&o,MapEdges);
                Quadrants.push_back(o);
            }
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * generateGridMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::setInitialState(vector<MeshPoint> &epts, vector<Quadrant> &eocts,
                                 map<QuadEdge, EdgeInfo> &medgs) {
        Quadrants.assign(make_move_iterator(eocts.begin()),
                         make_move_iterator(eocts.end()));
        
        //Erase previous quadrants to save memory
        eocts.erase(eocts.begin(),eocts.end());
        
        points.assign(make_move_iterator(epts.begin()),make_move_iterator(epts.end()));
        epts.erase(epts.begin(),epts.end());
        
        MapEdges = medgs;
        //recompute Octant indexes and update edges with those indexes
        for (unsigned int i=0;i<Quadrants.size();i++) {
            vector<unsigned int> qpts = Quadrants[i].getPointIndex();
            for (unsigned int j=0; j<4; j++) {
                
                auto qe = MapEdges.find(QuadEdge (qpts[j],qpts[(j+1)%4]));
                
                if (qe==MapEdges.end()) {
                    cout << "  Error at Mesher::setInitialState edge not found\n";
                }
                
                unsigned int mide = (qe->second)[0];
                if (mide!=0) {
                    auto qes1 = MapEdges.find(QuadEdge (qpts[j],mide));
                    auto qes2 = MapEdges.find(QuadEdge (qpts[(j+1)%4],mide));
                    if (j==0 || j==1) {
                        (qes1->second)[1] = i;
                        (qes2->second)[1] = i;
                    }
                    else {
                        (qes1->second)[2] = i;
                        (qes2->second)[2] = i;
                    }
                }
                else {
                    if (j==0 || j==1) {
                        (qe->second)[1] = i;
                    }
                    else {
                        (qe->second)[2] = i;
                    }
                }
            }
        }
        cout << "  Mesh initialized\n";

//        auto end_time = chrono::high_resolution_clock::now();
//        cout << "    * generateGridMesh in "
//        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
//        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::splitQuadrants(const unsigned short &rl, Polyline &input,
                                list<unsigned int> &roctli,
                                list<RefinementRegion *> &all_reg, const string &name,
                                const unsigned short &minrl,
                                const unsigned short &maxrl,
                                const bool &debugging){
        
        //The list of candidate quads to refine and the tmp version of
        //adding those how are still candidates for the next iteration.
        list<Quadrant> candidates, new_candidates, clean_processed, refine_tmp;
        
        //The Quads that don't need further refinement.
        vector<Quadrant> processed;
        
        //A map that connect the index of processed Quads with their position in the
        //vector of processed Quads.
        map<unsigned int, unsigned int> idx_pos_map;
        
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        
        //list<Quadrant>::iterator iter;
        
        //A set containing the index of Quads to be refined to
        //maintain balancing
        list<pair<unsigned int,unsigned int> > toBalance;
        
        //initialising the vector and the map
        candidates.assign(make_move_iterator(Quadrants.begin()),
                          make_move_iterator(Quadrants.end()));
        
        //The starting point for assignaiting indexes
        unsigned int new_q_idx = candidates.size();
        unsigned int future_q_idx = new_q_idx + 4*(unsigned int)roctli.size();
        
        //Erase previous quadrants to save memory
        Quadrants.erase(Quadrants.begin(),Quadrants.end());
        
        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(toBalance);
        
        auto start_refine_quad_time = chrono::high_resolution_clock::now();
        
        bool listref = false, localmax = false;
        if (!roctli.empty()) {
            listref = true;
        }
        
        
        //----------------------------------------------------------
        //refine once each Quadrant in the list
        //----------------------------------------------------------
        
        //note: roctli is an index list sorted when readed.
        //First we need to build the map with actual index
        //of quads in the processed quad vector. The quads
        //needing refinement will go into another list.
        if (!roctli.empty()) {
            
            list<unsigned int>::iterator octidx = roctli.begin();
            
            unsigned int qua_pos = 0;
            while (!candidates.empty()) {
                
                Quadrant quad = *(candidates.begin());
                candidates.pop_front();
                
                if (octidx == roctli.end() || qua_pos!=*octidx) {
                    idx_pos_map[quad.getIndex()] = processed.size();
                    processed.push_back(quad);
                }
                else {
                    //we advance to next quadrant in the list for the next iteration.
                    octidx++;
                    refine_tmp.push_back(quad);
                }
                qua_pos++;
            }
            
            for (auto quad: refine_tmp) {
                
                //start refinement process for current quadrant.
                list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                unsigned short qrl = quad.getRefinementLevel();
                
                if (qrl==maxrl) {
                    localmax = true;
                }
                
                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);
                
                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);
                sv.setStartIndex(new_q_idx);
                
                quad.accept(&sv);
                
                if (inter_edges.empty()) {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                        new_candidates.push_back(o);
                    }
                }
                else {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                        //the new points are inserted in bash at the end of this
                        //iteration. For this reason, the coordinates must be passed
                        //"manually" at this point (clipping_coords).
                        
                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_edges);
                        iv.setCoords(clipping_coords[j]);
                        
                        if (o.accept(&iv)) {
                            new_candidates.push_back(o);
                        }
                        else {
                            //The element doesn't intersect any input edge.
                            //It must be checked if it's inside or outside.
                            //Only in the first case add it to new_Quadrants.
                            //Test this with parent Quadrant faces only.
                            if (isItIn(input,inter_edges,clipping_coords[j])) {
                                new_candidates.push_back(o);
                            }
                            else {
                                //we must update neighbor information at the edges
                                EdgeVisitor::removeEdges(&o, MapEdges);
                            }
                        }
                    }
                }
            }
        
            //Erase the list to refine
            refine_tmp.erase(refine_tmp.begin(),refine_tmp.end());
            
            //cout << "To balance list has " << toBalance.size() << " quads\n";
        
            while (!toBalance.empty()) {
                
                list<pair<unsigned int, unsigned int> > tmp_toBalance;
                std::swap(toBalance,tmp_toBalance);
                tmp_toBalance.sort();
                tmp_toBalance.unique();
                
                while (!tmp_toBalance.empty()) {
                    unsigned int key = tmp_toBalance.begin()->first;
                    unsigned int val = tmp_toBalance.begin()->second;
                    
                    Quadrant quad = processed[val];
                    tmp_toBalance.pop_front();
                    list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                    unsigned short qrl = quad.getRefinementLevel();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    sv.setStartIndex(new_q_idx);
                    
                    quad.accept(&sv);
                    
                    if (inter_edges.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                idx_pos_map[o.getIndex()] = processed.size();
                                processed.push_back(o);
                                
                            }
                            else {
                                //The element doesn't intersect any input edge.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    idx_pos_map[o.getIndex()] = processed.size();
                                    processed.push_back(o);
                                }
                                else {
                                    //we must update neighbor information at the edges
                                    EdgeVisitor::removeEdges(&o, MapEdges);
                                }
                            }
                        }
                    }
                    //To mantain congruency in the map, we must erase all
                    //Quadrants (index) that have been split due to balancing.
                    auto delquad = idx_pos_map.find(key);
                    idx_pos_map.erase(delquad);
                }
            }
            
            // don't forget to update list
            std::swap(candidates,new_candidates);
            
            if (!new_pts.empty()) {
                //add the new points to the vector
                points.reserve(points.size() + new_pts.size());
                points.insert(points.end(),new_pts.begin(),new_pts.end());
                
                //cout << "new points inserted\n";
                //cout << "processed size " << processed.size() << endl;
                //cout << "map size " << idx_pos_map.size() << endl;
            }
            
            //unsigned int pro_quads = 0;
            
            //clean non used Quads.
            for (auto used_quad: idx_pos_map) {
                clean_processed.push_back(processed[used_quad.second]);
                //pro_quads++;
            }

            processed.erase(processed.begin(),processed.end());
            
        }
        
        //If there are no more refinement regions, we must apply transition patterns
        //at this moment and then finish the process.

        //----------------------------------------------------------
        // apply transition patterns
        //----------------------------------------------------------
        
        /*auto end_refine_quad_time = chrono::high_resolution_clock::now();
        
        //TransitionPatternVisitor section
        TransitionPatternVisitor tpv;
        tpv.setMapEdges(MapEdges);
        if (localmax) {
            tpv.setMaxRefLevel(maxrl+1);
        }
        else {
            tpv.setMaxRefLevel(maxrl);
        }
        new_pts.clear();

        //Apply transition patterns to remaining Quads
        unsigned mixedn = 0;
        for (auto &tq: clean_processed) {
            unsigned int sen = tq.getSubElements().size();
            if (!tq.accept(&tpv)) {
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
            if (tq.getSubElements().size()!=sen) {
                mixedn++;
            }
        }
        
        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }*/
        
        //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(candidates.begin()),make_move_iterator(candidates.end()));
        // better to erase as let in a indeterminate state by move
        candidates.erase(candidates.begin(),candidates.end());
        Quadrants.insert(Quadrants.end(),make_move_iterator(clean_processed.begin()),make_move_iterator(clean_processed.end()));
        clean_processed.erase(clean_processed.begin(),clean_processed.end());
        
#if (VTKOUT==true)            //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> transition_octree=make_shared<FEMesh>();
            saveOutputMesh(transition_octree,points,Quadrants);
            string tmp_name = name + "_transition";
            Services::WriteVTK(tmp_name,transition_octree);
        }
#endif
        
        if (localmax) {
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl+1,debugging,future_q_idx);
        }
        else {
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl,debugging,future_q_idx);
        }
        
        //If there are more refinement regions, continue with the process
        //were the quads positions will change in the final vector due to
        //map indexing that allows to optimize the research of neighbors.
        /*if (all_reg.size()>1) {
            
            auto end_time = chrono::high_resolution_clock::now();
            cout << "       * List refinement in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_refine_quad_time).count();
            cout << " ms"<< endl;
            
            if (listref) {
                
                //CL Debbuging
                {
                    //save pure octree mesh
                    std::shared_ptr<FEMesh> refined_octree=make_shared<FEMesh>();
                    saveOutputMesh(refined_octree,points,Quadrants);
                    string tmp_name = name + "_listRefinement";
                    Services::WriteVTK(tmp_name,refined_octree);
                }
            }
            //Continue with the rest of the refinement and apply transition patterns
            generateQuadtreeMesh(rl,input,all_reg,name,minrl,maxrl,debugging);
            return;
        }
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_refine_quad_time).count();
        cout << " ms"<< endl;
        cout << "    * generateQuadtreeMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_refine_quad_time).count();
        cout << " ms"<< endl;*/
        
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::generateQuadtreeMesh(const unsigned short &rl, Polyline &input,
                                      const list<RefinementRegion *> &all_reg,
                                      const string &name, const unsigned short &minrl,
                                      const unsigned short &givenmaxrl,
                                      const bool &debugging,
                                      unsigned int new_q_idx,bool Aliasing) {
        
        
        auto start_time = chrono::high_resolution_clock::now();
        
        //The list of candidate quads to refine and the tmp version of
        //adding those how are still candidates for the next iteration.
        list<Quadrant> candidates, new_candidates, refine_tmp, clean_processed;
        
        //The Quads that don't need further refinement.
        vector<Quadrant> processed;
        
        //A map that connect the index of processed Quads with their position in the
        //vector of processed Quads. 
        map<unsigned int, unsigned int> idx_pos_map;
        
        //list of the points added at this refinement iteration:
        list<Point3D> new_pts;
        
        //A set containing the index of Quads to be refined to
        //maintain balancing
        list<pair<unsigned int,unsigned int> > toBalance;
        
        //initialising the vector and the map
        candidates.assign(make_move_iterator(Quadrants.begin()),
                          make_move_iterator(Quadrants.end()));
        
        //Erase previous quadrants to save memory
        Quadrants.clear();
        
        //an iterator for the regions
        list<RefinementRegion *>::const_iterator reg_iter;
        
        //create visitors and give them variables
        SplitVisitor sv;
        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(toBalance);
        
        auto start_refine_quad_time = chrono::high_resolution_clock::now();
        
        //------------------------------------------------------------
        //refine each Quadrant until the Boundary is correctly handled
        //------------------------------------------------------------
        
        auto start_refine_rl_time = chrono::high_resolution_clock::now();
        
        //        list<RefinementRegion *>::const_iterator reg_iter=all_reg.begin();;
        unsigned int i=0; //current quad level
        /*do { // until no new quads are created
            new_pts.clear();
            
            //split the Quadrants as needed
            while (!tmp_Quadrants.empty()) {
                iter=tmp_Quadrants.begin();
                
                bool to_refine = false;
                
                iter->computeMaxDistance(points); //TODO, avoid recompute if already checked
                if ((*all_reg.begin())->intersectsQuadrant(points,*iter)) {
                    to_refine = true;
                }
                
                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    new_Quadrants.push_back(*iter);
                    tmp_Quadrants.pop_front();
                    continue;
                }
                else {
                    list<unsigned int> &inter_edges = iter->getIntersectedEdges();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    
                    sv.setStartIndex(new_Quadrants.size());
                    
                    //iter->split(points,new_pts,QuadEdges,split_elements,clipping_coords);
                    //cout << "Accept" << endl;
                    iter->accept(&sv);
                    
                    if (inter_edges.empty()) { //inner quad
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            new_Quadrants.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],i+1);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                new_Quadrants.push_back(o);
                            }
                            else {
                                //The element doesn't intersect any input face.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                
                                //Comment the following lines of this 'else' if
                                //only intersecting Quadrants are meant to be
                                //displayed.
                                
                                //note: inter_edges is quite enough to check if
                                //element is inside input, no Quadrant needed,
                                //so i moved the method to mesher  --setriva
                                
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    new_Quadrants.push_back(o);
                                }
                            }
                        }
                    }
                }
                // remove yet processed Quad
                tmp_Quadrants.pop_front();
                
            } // while
            
            // don't forget to update list
            std::swap(tmp_Quadrants,new_Quadrants);
            
            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cout << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                break;
            }
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
            
            
            ++i;
        } while (!new_pts.empty());
        
        #if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> bound_octree=make_shared<FEMesh>();
            saveOutputMesh(bound_octree,points,new_candidates);
            string tmp_name = name + "_bound";
            Services::WriteVTK(tmp_name,bound_octree);
        }
        #endif */

        unsigned short max_rl;
        if (givenmaxrl==0) {
            max_rl = i;
        }
        else {
            max_rl = givenmaxrl;
        }
        
        //cout << "rl : [" << minrl << "," << max_rl << "]\n";
        
        
        auto end_refine_rl_time = chrono::high_resolution_clock::now();
        cout << "         * boundary max " << i << " in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
        cout << " ms"<< endl;
        
        //----------------------------------------------------------
        //refine each Quadrant until the Refinement Level is reached
        //----------------------------------------------------------
        
        //when producing a new mesh, minrl will be always 0. But as this
        //method can also be called from "refineMesh" in this last case
        //the starting refinement level will not be 0, but the min rl
        //among the quadrants in the starting mesh.
        for (unsigned short i=minrl; i<rl; i++) {
            auto start_refine_rl_time = chrono::high_resolution_clock::now();
            
            //the new_pts is a list that holds the coordinates of
            //new points inserted at this iteration. At the end of
            //this bucle, they are inserted in the point vector
            new_pts.clear();
            
            //First we need to build the map with actual index
            //of quads in the processed quad vector. The quads
            //needing refinement will go into another list.
            while (!candidates.empty()) {
                Quadrant quad = *(candidates.begin());
                candidates.pop_front();
                
                bool to_refine = false;
                
                for (reg_iter=all_reg.begin(),reg_iter++; reg_iter!=all_reg.end(); ++reg_iter) {
                    
                    unsigned short region_rl = (*reg_iter)->getRefinementLevel();
                    if (region_rl<i) {
                        continue;
                    }
                    
                    //If the Quadrant has a greater RL than the region needs, continue
                    if (region_rl<= quad.getRefinementLevel()) {
                        continue;
                    }
                    
                    //in the case of a Region Refinement (a surface) this will change
                    //the state of inregion = true.
                    if ((*reg_iter)->intersectsQuadrant(points,quad)) {
                        to_refine = true;
                    }
                }
                
                //now if refinement is not needed, we add the Quadrant as it was.
                if (!to_refine) {
                    idx_pos_map[quad.getIndex()] = processed.size();
                    processed.push_back(quad);
                    continue;
                }
                else {
                    refine_tmp.push_back(quad);
                }
            }
            
            //now we can start to refine those needing it.
            for (auto quad: refine_tmp) {
                list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                unsigned short qrl = quad.getRefinementLevel();
                
                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);
                
                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);
                sv.setStartIndex(new_q_idx);
                
                quad.accept(&sv);
                
                //If quad has not intersected faces directly add sons to new
                //candidates
                if (inter_edges.empty()) {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                        //o.setInRegionState(reg_state);
                        new_candidates.push_back(o);
                    }
                }
                else {
                    for (unsigned int j=0; j<split_elements.size(); j++) {
                        Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                        //the new points are inserted in bash at the end of this
                        //iteration. For this reason, the coordinates must be passed
                        //"manually" at this point (clipping_coords).
                        
                        //select_faces = true
                        IntersectionsVisitor iv(true);
                        //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                        iv.setPolyline(input);
                        iv.setEdges(inter_edges);
                        iv.setCoords(clipping_coords[j]);
                        
                        if (o.accept(&iv)) {
                            //o.setInRegionState(reg_state);
                            new_candidates.push_back(o);
                        }
                        else {
                            //The element doesn't intersect any input edge.
                            //It must be checked if it's inside or outside.
                            //Only in the first case add it to new_Quadrants.
                            //Test this with parent Quadrant faces only.
                            if (isItIn(input,inter_edges,clipping_coords[j])) {
                                //o.setInRegionState(reg_state);
                                new_candidates.push_back(o);
                            }
                            else {
                                //we must update neighbor information at the edges
                                EdgeVisitor::removeEdges(&o, MapEdges);
                            }
                        }
                    }
                }
            }
            
            auto end_refine_rl_time = chrono::high_resolution_clock::now();
            cout << "         * level " << i << " in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_refine_rl_time-start_refine_rl_time).count();
            cout << " ms";
            
            //Erase the list to refine
            refine_tmp.erase(refine_tmp.begin(),refine_tmp.end());
        
            //Refine non balanced Quads
            while (!toBalance.empty()) {
                
                list<pair<unsigned int, unsigned int> > tmp_toBalance;
                std::swap(toBalance,tmp_toBalance);
                tmp_toBalance.sort();
                tmp_toBalance.unique();
                
                while (!tmp_toBalance.empty()) {
                    unsigned int key = tmp_toBalance.begin()->first;
                    unsigned int val = tmp_toBalance.begin()->second;
                    tmp_toBalance.pop_front();

                    auto delquad = idx_pos_map.find(key);
                    if (delquad == idx_pos_map.end()) {
                        //Note: let us say that quad N is in the tmp_to_balance list at current
                        //iteration. Before arriving to it, a neighbor of N is refined to such
                        //level that adds N to list to_balance (for the next iteration). However
                        //as N is in the current tmp_to_balance it will be removed from the
                        //map and in next iteration it will be attempt to be split and removed
                        //once more in the mesh. Therefore, this if avoids this case.
                        continue;
                    }
                    
                    Quadrant quad = processed[val];
                    list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                    unsigned short qrl = quad.getRefinementLevel();
                    
                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);
                    
                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    sv.setStartIndex(new_q_idx);
                    
                    quad.accept(&sv);
                    
                    if (inter_edges.empty()) {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j], qrl+1, new_q_idx++);
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    }
                    else {
                        for (unsigned int j=0; j<split_elements.size(); j++) {
                            Quadrant o (split_elements[j],qrl+1,new_q_idx++);
                            //the new points are inserted in bash at the end of this
                            //iteration. For this reason, the coordinates must be passed
                            //"manually" at this point (clipping_coords).
                            
                            //select_faces = true
                            IntersectionsVisitor iv(true);
                            //if (o.checkIntersections(input,inter_edges,clipping_coords[j]))
                            iv.setPolyline(input);
                            iv.setEdges(inter_edges);
                            iv.setCoords(clipping_coords[j]);
                            
                            if (o.accept(&iv)) {
                                idx_pos_map[o.getIndex()] = processed.size();
                                processed.push_back(o);
                                
                            }
                            else {
                                //The element doesn't intersect any input edge.
                                //It must be checked if it's inside or outside.
                                //Only in the first case add it to new_Quadrants.
                                //Test this with parent Quadrant faces only.
                                if (isItIn(input,inter_edges,clipping_coords[j])) {
                                    idx_pos_map[o.getIndex()] = processed.size();
                                    processed.push_back(o);
                                }
                                else {
                                    //we must update neighbor information at the edges
                                    EdgeVisitor::removeEdges(&o, MapEdges);
                                }
                            }
                        }
                    }
                    //To mantain congruency in the map, we must erase all
                    //Quadrants (index) that have been split due to balancing.
                    if (delquad != idx_pos_map.end()) {
                        idx_pos_map.erase(delquad);
                    }
                }
            }
            
            // don't forget to update list
            std::swap(candidates,new_candidates);
            
            //if no points were added at this iteration, it is no longer
            //necessary to continue the refinement.
            if (new_pts.empty()) {
                cerr << "warning at Mesher::generateQuadtreeMesh no new points!!!\n";
                auto end_balanced_time = chrono::high_resolution_clock::now();
                cout << " - Balanced in "
                << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_rl_time).count();
                cout << " ms"<< endl;
                break;
            }
            
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
            
            auto end_balanced_time = chrono::high_resolution_clock::now();
            cout << " - Balanced in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_balanced_time-end_refine_rl_time).count();
            cout << " ms"<< endl;
        }
        
        //clean non used Quads.
        for (auto used_quad: idx_pos_map) {
            clean_processed.push_back(processed[used_quad.second]);
        }
        
        //save space erasing processed quads.
        processed.erase(processed.begin(),processed.end());
        
#if (VTKOUT==true)        //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> refined_octree=make_shared<FEMesh>();
            if (debugging) {
                saveOutputMesh(refined_octree,points,clean_processed,true);
            }
            else {
                saveOutputMesh(refined_octree,points,clean_processed);
            }
            string tmp_name = name + "_refined";
            Services::WriteVTK(tmp_name,refined_octree);
        }
#endif

        auto end_refine_quad_time = chrono::high_resolution_clock::now();
        cout << "       * Refine Quad in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_refine_quad_time-start_refine_quad_time).count();
        cout << " ms"<< endl;

        
        //----------------------------------------------------------
        // apply transition patterns
        //----------------------------------------------------------
        /**/
        //TransitionPatternVisitor section
        TransitionPatternVisitor tpv;
        tpv.setMapEdges(MapEdges);
        tpv.setMaxRefLevel(max_rl);
        new_pts.clear();
        
        //Apply transition patterns to remaining Quads
        unsigned mixedn = 0;
        for (auto &tq: clean_processed) {
            unsigned int sen = tq.getSubElements().size();
            if (!tq.accept(&tpv)) {
                std::cerr << "Error at Mesher::generateQuadtreeMesh";
                std::cerr << " Transition Pattern not found\n";
            }
            if (tq.getSubElements().size()!=sen) {
                mixedn++;
            }
        }
        
        //if no points were added at this iteration, it is no longer
        //necessary to continue the refinement.
        if (!new_pts.empty()) {
            //add the new points to the vector
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(),new_pts.begin(),new_pts.end());
        }
        
        //insert will reserve space as well
        Quadrants.insert(Quadrants.end(),make_move_iterator(candidates.begin()),make_move_iterator(candidates.end()));
        // better to erase as let in a indeterminate state by move
        candidates.erase(candidates.begin(),candidates.end());
        Quadrants.insert(Quadrants.end(),make_move_iterator(clean_processed.begin()),make_move_iterator(clean_processed.end()));
        clean_processed.erase(clean_processed.begin(),clean_processed.end());
        
#if (VTKOUT==true) //CL Debbuging
        {
            //save pure octree mesh
            std::shared_ptr<FEMesh> transition_octree=make_shared<FEMesh>();
            saveOutputMesh(transition_octree,points,Quadrants);
            string tmp_name = name + "_transition";
            Services::WriteVTK(tmp_name,transition_octree);
        }
#endif

        auto end_time = chrono::high_resolution_clock::now();
        cout << "       * Transition Patterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-end_refine_quad_time).count();
        cout << " ms"<< endl;
        cout << "    * generateQuadtreeMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    bool Mesher::isItIn(const Polyline &mesh, const list<unsigned int> &faces, const vector<Point3D> &coords) const {
        //this method is meant to be used by Quadrants that don't
        //intersect input domains. If they are inside of at least
        //one input mesh, then they must remain in the output mesh.
        
        bool first = mesh.pointIsInMesh(coords[0],faces);
        bool second = mesh.pointIsInMesh(coords[1],faces);
        if (first==second) {
            return first;
        }
        
        //cout << "one inconsistency detected -> hard test\n";
        //return mesh.pointIsInMesh(coords[0],faces);
        return mesh.pointIsInMesh(coords[0]);
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh, bool decoration){
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<vector<unsigned int> > out_els;
        vector<unsigned short > out_els_ref_level, out_els_surf, out_els_deb;
        vector<double > out_els_min_angle, out_els_max_angle;
        array <unsigned int,180> out_els_angle_tri_histogram = {0},
                out_els_angle_quad_histogram = {0};

        //even if we don't know the quantity of elements, it will be at least the
        //number of quadrants, so we reserve for the vectors of VTK decoration.
        out_els_min_angle.reserve(Quadrants.size());
        out_els_max_angle.reserve(Quadrants.size());
        out_els_ref_level.reserve(Quadrants.size());
        out_els_surf.reserve(Quadrants.size());
        
        //new_idxs will hold the index of used nodes in the outside vector for points.
        //If the a node is not used by any element, its index will be 0 in this vector,
        //therefore the actual index is shiffted in 1. In other words, node 0 is node 1,
        //and node n is node n+1.
        vector<unsigned int> new_idxs (points.size(),0);
        unsigned int out_node_count = 0;
        vector<Point3D> out_pts;
 
        //recompute node indexes and update elements with them.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            const vector<vector<unsigned int> > &sub_els= Quadrants[i].getSubElements();
            for (unsigned int j=0; j<sub_els.size(); j++) {
                
                vector<unsigned int> sub_ele_new_idxs = sub_els[j];
                for (unsigned int k=0; k<sub_ele_new_idxs.size();k++) {
                    
                    unsigned int p_idx = sub_ele_new_idxs[k];
                    
                    if (new_idxs[p_idx]==0) {
                        sub_ele_new_idxs[k] = out_node_count++;
                        new_idxs[p_idx]=out_node_count;
                        out_pts.push_back(points[p_idx].getPoint());
                    }
                    else {
                        sub_ele_new_idxs[k] = new_idxs[p_idx]-1;
                    }
                }
                if (decoration) {
                    //surface regarding quad
                    if (Quadrants[i].isSurface()) {
                        out_els_surf.push_back(1);
                    }
                    else {
                        out_els_surf.push_back(0);
                    }
                    
                    if (Quadrants[i].isDebugging()) {
                        out_els_deb.push_back(1);
                    }
                    else {
                        out_els_deb.push_back(0);
                    }
                    
                    //refinment level herited from quad
                    out_els_ref_level.push_back(Quadrants[i].getRefinementLevel());


                    //compute minAngle and maxAngle (only for elements on the surface)
                    unsigned int np=sub_ele_new_idxs.size(); //nb points of the element
                    double minAngle=std::numeric_limits<double>::max();
                    double maxAngle=std::numeric_limits<double>::min();
                    if (Quadrants[i].isSurface()) {
                        for (unsigned int k=0; k<np; ++k) {

                            const Point3D &P0 = out_pts[sub_ele_new_idxs[(k-1+np)%np]];
                            const Point3D &P1 = out_pts[sub_ele_new_idxs[k]];
                            const Point3D &P2 = out_pts[sub_ele_new_idxs[(k+1)%np]];

                            double angle=P1.angle3Points(P0,P2);
                            minAngle=std::min(minAngle, angle);
                            maxAngle=std::max(maxAngle, angle);
                            if (np==3)
                                out_els_angle_tri_histogram[ min(179,((int)( round(angle)/*/10.*/)) %180) ]++;
                            else if (angle<89.5 || angle>91.5)
                                out_els_angle_quad_histogram[ min(179,((int)(round(angle)/*/10.*/)) %180) ]++;
                        }
                    }
                    out_els_min_angle.push_back(minAngle);
                    out_els_max_angle.push_back(maxAngle);
                }
                out_els.push_back(sub_ele_new_idxs);
            }
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        mesh->setRefLevels(out_els_ref_level);
        mesh->setMinAngles(out_els_min_angle);
        mesh->setMaxAngles(out_els_max_angle);
        mesh->setAnglesTriHistogram(out_els_angle_tri_histogram);
        mesh->setAnglesQuadHistogram(out_els_angle_quad_histogram);
        mesh->setSurfState(out_els_surf);
        mesh->setDebugging(out_els_deb);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                        vector<MeshPoint> &tmp_points,
                                        list<Quadrant> &tmp_Quadrants,
                                        const bool &debugging,
                                        const list<Point3D> &extra_pts) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;
        vector<unsigned short> deb_els;
        
        if (debugging) {
            deb_els.reserve(tmp_Quadrants.size());
        }
        
        unsigned int n = tmp_points.size();
        out_pts.reserve(n + extra_pts.size());
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }
        
        //copy tmp points if available
        for (auto p: extra_pts) {
            out_pts.push_back(p);
        }
        
        OneIrregularVisitor oiv;
        if (debugging) {
            oiv.setEdges(MapEdges);
        }
        
        list<Quadrant>::iterator o_iter;
        
        bool errors = false;
        if (debugging) {
            for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {
                
                vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
                //Check balancing
                bool result = o_iter->accept(&oiv);
                if (result) {
                    deb_els.push_back(0);
                }
                else {
                    errors = true;
                    deb_els.push_back(1);
                }
            }
            mesh->setDebugging(deb_els);
            if (errors) {
                cerr << "->Warning mesh not balanced at Mesher::saveOutputMesh";
                cerr << " with debugging option on.\n";
            }
            else {
                cout << "    * No balancing problems detected\n";
            }
        }
        else {
            for (o_iter=tmp_Quadrants.begin(); o_iter!=tmp_Quadrants.end(); ++o_iter) {
                
                vector<vector<unsigned int> > sub_els= o_iter->getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
            }
        }
        
        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;
        
        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    unsigned int Mesher::saveOutputMesh(const shared_ptr<FEMesh> &mesh,
                                        vector<MeshPoint> &tmp_points,
                                        vector<Quadrant> &tmp_Quadrants,
                                        const bool &debugging,
                                        const list<Point3D> &extra_pts) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        vector<Point3D> out_pts;
        list<vector<unsigned int> > tmp_elements;
        vector<vector<unsigned int> > out_els;
        vector<unsigned short> deb_els;
        
        if (debugging) {
            deb_els.reserve(tmp_Quadrants.size());
        }
        
        unsigned int n = tmp_points.size();
        out_pts.reserve(n+extra_pts.size());
        for (unsigned int i=0; i<n; i++) {
            out_pts.push_back(points[i].getPoint());
        }
        
        OneIrregularVisitor oiv;
        if (debugging) {
            oiv.setEdges(MapEdges);
        }
        
        
        //copy tmp points if available
        for (auto p: extra_pts) {
            out_pts.push_back(p);
        }

        // Guard: if all quadrants were dropped (e.g. by the
        // archipelago resolver at L >= all component sizes, or by a
        // future caller that empties Quadrants), emit an empty mesh
        // rather than dereferencing tmp_Quadrants[0]. This was a
        // pre-existing crash reachable from the post-resolve
        // debugging/projection/removal paths at lines 138/150/165.
        if (tmp_Quadrants.empty()) {
            mesh->setPoints(out_pts);
            mesh->setElements(out_els);
            auto end_time = chrono::high_resolution_clock::now();
            cout << "    * SaveOutputMesh (empty Quadrants) in "
                 << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count()
                 << " ms" << endl;
            return 0;
        }

        unsigned int frl = tmp_Quadrants[0].getRefinementLevel();

        bool errors = false;
        
        if (debugging) {
            for (unsigned int i=0; i<tmp_Quadrants.size(); i++) {
                vector<vector<unsigned int> > sub_els= tmp_Quadrants[i].getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
                //Check balancing
                bool result = tmp_Quadrants[i].accept(&oiv);
                if (result) {
                    deb_els.push_back(0);
                }
                else {
                    errors = true;
                    deb_els.push_back(1);
                }
            }
            mesh->setDebugging(deb_els);
            if (errors) {
                cerr << "->Warning mesh not balanced at Mesher::saveOutputMesh";
                cerr << " with debugging option on.\n";
            }
            else {
                cout << "    * No balancing problems detected\n";
            }
        }
        else {
            for (unsigned int i=0; i<tmp_Quadrants.size(); i++) {
                vector<vector<unsigned int> > sub_els= tmp_Quadrants[i].getSubElements();
                for (unsigned int j=0; j<sub_els.size(); j++) {
                    tmp_elements.push_back(sub_els[j]);
                }
            }
        }
        
        out_els.reserve(tmp_elements.size());
        list<vector<unsigned int> >::const_iterator e_iter;
        
        for (e_iter=tmp_elements.begin(); e_iter!=tmp_elements.end(); ++e_iter) {
            out_els.push_back(*e_iter);
        }
        
        mesh->setPoints(out_pts);
        mesh->setElements(out_els);
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * SaveOutputMesh in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
        return out_els.size();
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::detectFeatureQuadrants(Polyline &input) {
        
        auto start_time = chrono::high_resolution_clock::now();
        
        unsigned int featCount = 0;
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (input.getNbFeatures(Quadrants[i],points)>0) {
                //Quadrants[i].setFeature(); now set by Polyline::getNbFeatures()
                featCount++;
            }
        }
        //cout << "number of quadrants with features: " << featCount << "\n";
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectFeatureQuadrants in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
        
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::linkElementsToNodes(){
        
        auto start_time = chrono::high_resolution_clock::now();
        
        //clear previous information
        for (unsigned int i=0; i<points.size(); i++) {
            points[i].clearElements();
        }
        
        //link element info to nodes
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            
            const vector <unsigned int> &q_indpts = Quadrants[i].getSubPointIndex();
            
            for (unsigned int j=0; j<q_indpts.size(); j++) {
                points.at(q_indpts[j]).addElement(i);
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * linkElementsToNodes in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    //WARNING: linkElementsToNodes() must be called before
    void Mesher::computeNodeMaxDist() {
        for (auto &q:Quadrants) {
            q.computeMaxDistance(points);
        }
        for (auto& p:points) {
            for (const auto &e:p.getElements()) {
                p.updateMaxDistance(Quadrants[e].getMaxDistance());
            }
        }
}

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    // TUSQH-style subdivision loop.
    //
    // Each iteration:
    //   1. Run WindingNumberSubdivisionVisitor on every candidate quadrant.
    //      The visitor:
    //         - computes the s x s winding numbers;
    //         - classifies the quadrant as AllInside / AllOutside / Mixed;
    //         - returns true iff the quadrant must be subdivided.
    //   2. Quadrants classified as AllInside or AllOutside are pushed to
    //      `processed` (no further refinement).
    //   3. Quadrants classified as Mixed (or forced by the optional
    //      edge-intersect criterion) are split via SplitVisitor. Their 4
    //      children become the new `candidates` for the next iteration.
    //   4. After splitting, balance the quadtree as usual using toBalance /
    //      idx_pos_map.
    //
    // Stops when:
    //   - maxDepth is reached, or
    //   - no candidate was subdivided in this iteration.
    //
    // This mirrors the "subdivide while ambiguous" rule from the TUSQH
    // paper.
    void Mesher::windingSubdivide(Polyline &input, unsigned int maxDepth,
                                  unsigned int tusqhSampleSize,
                                  bool refineOnEdgeIntersect,
                                  const string &name,
                                  const bool &debugging,
                                  unsigned int tusqhExtraResolveDepth,bool Aliasing)
    {
        auto start_time = chrono::high_resolution_clock::now();
        
        // The list of candidate quadrants to refine and the temporary
        // list of new candidates created at the current iteration.
        list<Quadrant> candidates, new_candidates;

        // Quadrants that don't need further refinement this iteration.
        vector<Quadrant> processed;

        // Map: quadrant index -> position in `processed` (used for
        // balancing lookups, like the rest of the Mesher code).
        map<unsigned int, unsigned int> idx_pos_map;

        // Points added at this iteration (mid-edge + center nodes).
        list<Point3D> new_pts;

        // Quadrant indices whose WindingState must be recomputed AFTER
        // this iteration's new_pts have been appended to `points`
        // (otherwise getSamplePoint() reads past the end of `points`).
        // Populated by the balance path; cleared after classification.
        vector<unsigned int> needs_classification;

        // Quad indices that need to be split to keep the quadtree
        // one-irregular.
        list<pair<unsigned int, unsigned int> > toBalance;

        // Seed: the root quadrants produced by generateGridMesh /
        // generateQuadtreeMesh (i.e. the cells that survived the
        // classical pipeline and are now the input to TUSQH). Their
        // origin is tagged as Classical so the post-TUSQH snapshots
        // can distinguish them from cells born inside the TUSQH loop.
        candidates.assign(make_move_iterator(Quadrants.begin()),
                          make_move_iterator(Quadrants.end()));
        Quadrants.clear();
        for (auto& q : candidates) {
            q.setOrigin(QuadrantOrigin::Classical);
        }

        // Visitors reused across iterations.
        WindingNumberSubdivisionVisitor wnsv(tusqhSampleSize,
                                              refineOnEdgeIntersect);
        wnsv.setPolyline(&input);
        wnsv.setPoints(&points);

        // We keep one IntersectionsVisitor around only when the legacy
        // criterion is requested.
        IntersectionsVisitor legacyIv(true);
        if (refineOnEdgeIntersect) {
            wnsv.setIntersectionsVisitor(&legacyIv);
        }

        SplitVisitor sv;
        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(toBalance);

        // Persistent index counter for new (child) quadrants. We never
        // reset this inside the depth loop, otherwise we'd produce
        // duplicate q_ids across iterations.
        unsigned int new_q_idx = 0;
        for (const auto& q : candidates) {
            unsigned int qi = q.getIndex();
            if (qi != std::numeric_limits<unsigned int>::max() && qi + 1 > new_q_idx) {
                new_q_idx = qi + 1;
            }
        }
        // Respaldo por si acaso:
        if (new_q_idx < candidates.size()) {
            new_q_idx = candidates.size();
        }

        // Calculamos el nivel de refinamiento (qrl) máximo de los candidatos
        // de entrada. Esto es necesario porque el parámetro `maxDepth`
        // representa el qrl MÁXIMO absoluto que deben alcanzar los hijos,
        // no el número de iteraciones del bucle. Si la grilla ya fue
        // pre-refinada (e.g. por preRefineForTusqh, que lleva la grilla
        // inicial a baseLevel) o si estamos refinando un .oct ya profundo
        // (caso refineMesh con useTusqh), debemos acotar el bucle para
        // que los hijos nunca excedan `maxDepth`.
        //
        // Casos cubiertos:
        //   startDepth == 0         -> effectiveMaxDepth == maxDepth
        //                              (caso simple, comportamiento idéntico al original)
        //   startDepth <  maxDepth  -> effectiveMaxDepth == maxDepth - startDepth
        //                              (caso Chesapeake Bay con heurística: arregla el bug)
        //   startDepth >= maxDepth  -> effectiveMaxDepth == 0
        //                              (bucle no corre; candidatos pasan al post-proceso)
        unsigned short startDepth = 0;
        for (const auto& q : candidates) {
            if (q.getRefinementLevel() > startDepth) {
                startDepth = q.getRefinementLevel();
            }
        }
        const unsigned int effectiveMaxDepth =
            (startDepth < maxDepth) ? (maxDepth - startDepth) : 0u;

        if (effectiveMaxDepth == 0 && !candidates.empty()) {
            cout << "    * windingSubdivide: candidatos ya en qrl="
                 << startDepth << " >= maxDepth=" << maxDepth
                 << "; el bucle TUSQH no subdividirá más.\n";
        }

        for (unsigned int depth = 0; depth < effectiveMaxDepth; ++depth) {
            auto start_depth_time = chrono::high_resolution_clock::now();
            new_pts.clear();
            new_candidates.clear();

            if (candidates.empty()) {
                break;
            }

            // 1) Classify every candidate using the TUSQH criterion.
            list<Quadrant> refine_tmp;
            refine_tmp.clear();

            unsigned int refinedCount = 0;

            // Index of the first entry in `processed` added by THIS
            // iteration. We use it to dump only the cells that were
            // classified this depth, leaving previous iterations'
            // entries alone (so the per-iter snapshot stays small and
            // meaningful).
            const size_t processed_size_before = processed.size();

            while (!candidates.empty()) {
                Quadrant quad = *(candidates.begin());
                candidates.pop_front();

                bool toRefine = wnsv.visit(&quad);

                if (!toRefine) {
                    // Inside or outside, no subdivision.
                    idx_pos_map[quad.getIndex()] = processed.size();
                    processed.push_back(quad);
                } else {
                    refine_tmp.push_back(quad);
                    refinedCount++;
                }
            }

#if (VTKOUT==true)
            {
                // Per-iteration debug dumps. After classification, every
                // candidate has been classified (WindingState set,
                // mWindingNumbers populated, mVolumeFraction computed)
                // and the to-be-split cells live in `refine_tmp`. We
                // dump two views of this iteration:
                //   - `<name>_tusqh_iter<N>.vtk`        quadtree snapshot
                //                                         of all classified
                //                                         cells with
                //                                         winding_state,
                //                                         VF, ref_level
                //                                         and sample_size.
                //   - `<name>_tusqh_iter<N>_samples.vtk` heat map of every
                //                                         s x s sample
                //                                         with its
                //                                         winding_number
                //                                         and parent
                //                                         q_id.
                // Together they let the user inspect in ParaView exactly
                // which cells were deemed Mixed (and why) at each depth,
                // and which samples pushed them over the threshold.
                string iter_label = "_tusqh_iter" + std::to_string(depth);
                string snap_name = name + iter_label;
                string samples_name = name + iter_label;

                // Collect this iteration's classified cells: the new
                // entries pushed to `processed` (AllInside / AllOutside
                // leaves) plus the entries in `refine_tmp` (Mixed cells
                // about to be split). Cells whose mWindingNumbers was
                // never populated (e.g. legacy `-E` IntersectionsVisitor
                // path) are skipped silently by the templated
                // writeVFHeatmap.
                std::vector<Quadrant> classified;
                classified.reserve((processed.size() - processed_size_before)
                                   + refine_tmp.size());
                for (size_t i = processed_size_before; i < processed.size(); ++i) {
                    classified.push_back(processed[i]);
                }
                for (const auto& q : refine_tmp) {
                    classified.push_back(q);
                }

                // Snapshot: <name>_tusqh_iter<N>.vtk
                VolumeFractionVTKWriter::writeCandidatesSnapshot(
                    snap_name, classified, points, "");
                // Heat map: <name>_tusqh_iter<N>_samples.vtk (the
                // function appends "_samples" to the filename).
                VolumeFractionVTKWriter::writeVFHeatmap(
                    samples_name, classified, points);
            }
#endif

            // 2) Stop if no candidate was subdivided (the quadtree is
            //    already "TUSQH-stable").
            if (refine_tmp.empty()) {
                auto end_depth_time = chrono::high_resolution_clock::now();
                cout << "         * TUSQH iter " << depth << "/" << effectiveMaxDepth
                     << " (target qrl=" << (startDepth + depth + 1)
                     << ") (no refinement needed) in "
                     << std::chrono::duration_cast<chrono::milliseconds>(end_depth_time-start_depth_time).count()
                     << " ms" << endl;
                // Early-stop path: do NOT clear processed here, because
                // the post-loop promotion below still expects to read
                // entries from processed via idx_pos_map.
                break;
            }

            // 3) Split the mixed cells.
            for (auto& quad : refine_tmp) {
                list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                unsigned short qrl = quad.getRefinementLevel();

                vector<vector<Point3D> > clipping_coords;
                sv.setClipping(clipping_coords);

                vector<vector<unsigned int> > split_elements;
                sv.setNewEles(split_elements);
                sv.setStartIndex(new_q_idx);

                quad.accept(&sv);

                if (inter_edges.empty()) {
                    // Pure interior split (no Polyline edge nearby).
                    for (unsigned int j = 0; j < split_elements.size(); j++) {
                        Quadrant o(split_elements[j], qrl+1, new_q_idx++);
                        o.setOrigin(QuadrantOrigin::TusqhSplit);
                        new_candidates.push_back(o);
                    }
                } else {
                    // TUSQH refinement of a boundary cell. Children are
                    // queued without re-running the IntersectionsVisitor:
                    // the next iteration will classify them again from
                    // scratch (their winding numbers may now be
                    // unambiguous).
                    for (unsigned int j = 0; j < split_elements.size(); j++) {
                        Quadrant o(split_elements[j], qrl+1, new_q_idx++);
                        o.setOrigin(QuadrantOrigin::TusqhSplit);
                        new_candidates.push_back(o);
                    }
                }
            }
            refine_tmp.clear();

            // 4) Balance the quadtree by splitting neighbours that
            //    became too coarse with respect to their refined
            //    siblings.
            while (!toBalance.empty()) {
                list<pair<unsigned int, unsigned int> > tmp_toBalance;
                std::swap(toBalance, tmp_toBalance);
                tmp_toBalance.sort();
                tmp_toBalance.unique();

                while (!tmp_toBalance.empty()) {
                    unsigned int key = tmp_toBalance.begin()->first;
                    unsigned int val = tmp_toBalance.begin()->second;
                    tmp_toBalance.pop_front();

                    auto delquad = idx_pos_map.find(key);
                    if (delquad == idx_pos_map.end()) {
                        continue;
                    }

                    Quadrant quad = processed[val];
                    list<unsigned int> &inter_edges = quad.getIntersectedEdges();
                    unsigned short qrl = quad.getRefinementLevel();

                    vector<vector<Point3D> > clipping_coords;
                    sv.setClipping(clipping_coords);

                    vector<vector<unsigned int> > split_elements;
                    sv.setNewEles(split_elements);
                    sv.setStartIndex(new_q_idx);

                    quad.accept(&sv);

                    if (inter_edges.empty()) {
                        for (unsigned int j = 0; j < split_elements.size(); j++) {
                            Quadrant o(split_elements[j], qrl+1, new_q_idx++);
                            o.setOrigin(QuadrantOrigin::TusqhBalance);
                            // Defer classification until points are
                            // appended to the global point list.
                            needs_classification.push_back(o.getIndex());
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    } else {
                        for (unsigned int j = 0; j < split_elements.size(); j++) {
                            Quadrant o(split_elements[j], qrl+1, new_q_idx++);
                            o.setOrigin(QuadrantOrigin::TusqhBalance);
                            needs_classification.push_back(o.getIndex());
                            idx_pos_map[o.getIndex()] = processed.size();
                            processed.push_back(o);
                        }
                    }

                    if (delquad != idx_pos_map.end()) {
                        idx_pos_map.erase(delquad);
                    }
                }
            }

            // Promote the refined children to be the next round's
            // candidates.
            std::swap(candidates, new_candidates);

            if (!new_pts.empty()) {
                points.reserve(points.size() + new_pts.size());
                points.insert(points.end(), new_pts.begin(), new_pts.end());
            }

            // Deferred classification: balanced children were pushed to
            // `processed` before their midpoint points existed. Now that
            // `points` includes the new nodes, classify them so the
            // AllOutside filter (below) can drop them.
            for (auto qid : needs_classification) {
                auto it = idx_pos_map.find(qid);
                if (it == idx_pos_map.end()) continue;
                wnsv.visit(&processed[it->second]);
            }
            needs_classification.clear();

            auto end_depth_time = chrono::high_resolution_clock::now();
            cout << "         * TUSQH iter " << depth << "/" << effectiveMaxDepth
                 << " (target qrl=" << (startDepth + depth + 1)
                 << ") (refined " << refinedCount << " cells) in "
                 << std::chrono::duration_cast<chrono::milliseconds>(end_depth_time-start_depth_time).count()
                 << " ms" << endl;
        }

        // Any candidate left over (e.g. when we hit maxDepth) is
        // promoted to processed, but we still need to classify them
        // (their WindingState may be Unknown from the constructor).
        for (auto& q : candidates) {
            wnsv.visit(&q);
            idx_pos_map[q.getIndex()] = processed.size();
            processed.push_back(q);
        }
        candidates.clear();

        // TUSQH post-processing: PRESERVE the full cubical complex.
        //
        // Paper-faithful (TUSQH §3.4): the bridge-joining algorithm
        // needs access to the background grid topology. Cells with
        // AllOutside winding state (0% volume fraction) are kept in
        // `Quadrants` so that MapEdges retains the edges between them
        // and the surrounding interior cells. This lets the bridge
        // algorithm identify pairs of components whose only separation
        // is a chain of exterior cells and connect them.
        //
        // The downstream pipeline treats AllOutside cells as invisible:
        //   - linkElementsToNodes skips them,
        //   - detectInsideNodes / detectFeatureQuadrants skip them,
        //   - removeOnSurfaceSafe already skips them via isInside(),
        //   - saveOutputMesh filters them out at output time.
        // The intersected_edges list of AllInside cells is cleared to
        // honour the classical pipeline contract
        // (isInside() == intersected_edges.empty()). Mixed and Unknown
        // cells are kept untouched with their inherited edge list.
        unsigned int preservedOutside = 0;
        for (auto& used_quad : idx_pos_map) {
            Quadrant& q = processed[used_quad.second];
            switch (q.getWindingState()) {
                case WindingState::AllOutside:
                    // Paper-faithful: keep the cell so the cubical
                    // complex remains intact for bridge-joining. Its
                    // isInside() returns false, so all downstream
                    // consumers skip it.
                    q.getIntersectedEdges().clear();
                    q.setOrigin(QuadrantOrigin::PreservedOutsider);
                    Quadrants.push_back(std::move(q));
                    ++preservedOutside;
                    break;
                case WindingState::AllInside:
                    q.getIntersectedEdges().clear();
                    Quadrants.push_back(std::move(q));
                    break;
                case WindingState::Mixed:
                case WindingState::Unknown:
                default:
                    Quadrants.push_back(std::move(q));
                    break;
            }
        }
        processed.clear();

        if (preservedOutside > 0) {
            cout << "    * TUSQH preserved " << preservedOutside
                 << " AllOutside cells in cubical complex (0% VF, "
                 << "used by bridge-joining for topology)\n";
        }

        // ----------------------------------------------------------------
        // TUSQH resolve pass (optional).
        //
        // After the main depth loop, some cells may still be Mixed at
        // qrl == maxDepth (they hit the user's hard refinement bound
        // before TUSQH could resolve them). If `tusqhExtraResolveDepth`
        // is positive, we run additional pure-TUSQH iterations on those
        // Mixed cells only, allowing their children to exceed maxDepth.
        //
        // Properties of this pass:
        //   - No balance: we only subdivide Mixed cells and accept the
        //     resulting one-irregular violations (the alternative would
        //     be a recursive resolve of neighbours, which is O(N^2)).
        //   - No legacy edge-intersect criterion: pure TUSQH criterion
        //     throughout.
        //   - AllOutside filter still applies: cells whose s x s samples
        //     are uniformly zero are dropped before reaching the output.
        //   - The depth bound is `tusqhExtraResolveDepth`, so the user
        //     stays in control of the maximum cell count.
        // ----------------------------------------------------------------
        if (tusqhExtraResolveDepth > 0) {
            unsigned int mixedBefore = 0;
            list<Quadrant> resolve_candidates;
            vector<Quadrant> kept_non_mixed;
            kept_non_mixed.reserve(Quadrants.size());
            for (auto& q : Quadrants) {
                if (q.getWindingState() == WindingState::Mixed) {
                    resolve_candidates.push_back(std::move(q));
                    ++mixedBefore;
                } else {
                    kept_non_mixed.push_back(std::move(q));
                }
            }
            Quadrants.clear();

            cout << "    * TUSQH resolve pass: " << mixedBefore
                 << " Mixed leaves at qrl=" << maxDepth << ", "
                 << tusqhExtraResolveDepth << " extra depth\n";

            if (!resolve_candidates.empty()) {
                auto resolve_start = chrono::high_resolution_clock::now();

                // Visitors reused across resolve iterations.
                WindingNumberSubdivisionVisitor resolveVisitor(tusqhSampleSize,
                                                               false);
                resolveVisitor.setPolyline(&input);
                resolveVisitor.setPoints(&points);

                list<Quadrant> new_resolve_candidates;
                vector<Quadrant> resolve_processed;
                list<Point3D> resolve_new_pts;
                vector<unsigned int> resolve_needs_classification;

                // Empty container references for SplitVisitor: we keep
                // it happy (its `proQuadMap` and `toBalanceList`
                // pointers need to be non-NULL because it queries them
                // unconditionally) but never push into them, so balance
                // is a no-op in the resolve pass.
                list<pair<unsigned int, unsigned int> > emptyToBalance;
                map<unsigned int, unsigned int> emptyIdxPosMap;

                SplitVisitor resolveSv;
                resolveSv.setPoints(points);
                resolveSv.setMapEdges(MapEdges);
                resolveSv.setNewPts(resolve_new_pts);
                resolveSv.setProcessedQuadVector(resolve_processed);
                resolveSv.setMapProcessed(emptyIdxPosMap);
                resolveSv.setToBalanceList(emptyToBalance);

                unsigned int totalResolved = 0;
                unsigned int totalDroppedInResolve = 0;

                for (unsigned int rdepth = 0; rdepth < tusqhExtraResolveDepth; ++rdepth) {
                    if (resolve_candidates.empty()) {
                        break;
                    }
                    auto dstart = chrono::high_resolution_clock::now();

                    resolve_new_pts.clear();
                    new_resolve_candidates.clear();
                    resolve_needs_classification.clear();
                    unsigned int refinedInIter = 0;

                    while (!resolve_candidates.empty()) {
                        Quadrant quad = *(resolve_candidates.begin());
                        resolve_candidates.pop_front();

                        bool toRefine = resolveVisitor.visit(&quad);
                        if (!toRefine) {
                            // AllInside or AllOutside: leaf. Will be
                            // filtered below.
                            resolve_processed.push_back(std::move(quad));
                        } else {
                            // Mixed: split into 4 children. No balance
                            // in the resolve pass.
                            list<unsigned int> &inter_edges =
                                quad.getIntersectedEdges();
                            unsigned short qrl = quad.getRefinementLevel();

                            vector<vector<Point3D> > clipping_coords;
                            resolveSv.setClipping(clipping_coords);
                            vector<vector<unsigned int> > split_elements;
                            resolveSv.setNewEles(split_elements);
                            resolveSv.setStartIndex(new_q_idx);

                            quad.accept(&resolveSv);

                            for (unsigned int j = 0; j < split_elements.size(); ++j) {
                                Quadrant o(split_elements[j], qrl+1, new_q_idx++);
                                o.setOrigin(QuadrantOrigin::TusqhResolve);
                                // Pure-TUSQH resolve: classify immediately
                                // using the s x s samples of the child.
                                // The child's mid-edge points were just
                                // appended by SplitVisitor into
                                // resolve_new_pts (not yet in `points`),
                                // so we defer classification to after
                                // resolve_new_pts is appended below.
                                resolve_needs_classification.push_back(o.getIndex());
                                new_resolve_candidates.push_back(std::move(o));
                            }
                            ++refinedInIter;
                        }
                    }

                    // Append mid-edge points created this iteration.
                    if (!resolve_new_pts.empty()) {
                        points.reserve(points.size() + resolve_new_pts.size());
                        points.insert(points.end(),
                                      resolve_new_pts.begin(),
                                      resolve_new_pts.end());
                    }

                    // Deferred classification of new children.
                    map<unsigned int, unsigned int> resolve_idx_pos_map;
                    for (auto& q : new_resolve_candidates) {
                        resolve_idx_pos_map[q.getIndex()] = resolve_processed.size();
                        resolve_processed.push_back(std::move(q));
                    }
                    for (auto qid : resolve_needs_classification) {
                        auto it = resolve_idx_pos_map.find(qid);
                        if (it == resolve_idx_pos_map.end()) continue;
                        resolveVisitor.visit(&resolve_processed[it->second]);
                    }
                    resolve_needs_classification.clear();

                    // Promote non-Mixed to processed; keep Mixed for next
                    // iteration.
                    list<Quadrant> next_iter_candidates;
                    for (auto& q : resolve_processed) {
                        if (q.getWindingState() == WindingState::Mixed) {
                            next_iter_candidates.push_back(std::move(q));
                        } else {
                            ++totalResolved;
                            if (q.getWindingState() == WindingState::AllOutside) {
                                ++totalDroppedInResolve;
                            }
                            // AllInside stays, its intersected_edges is
                            // already empty (no edge was assigned in the
                            // resolve pass because we did not run the
                            // IntersectionsVisitor here).
                        }
                    }
                    resolve_processed.clear();

                    std::swap(resolve_candidates, next_iter_candidates);

                    auto dend = chrono::high_resolution_clock::now();
                    cout << "      - resolve depth " << rdepth
                         << ": refined " << refinedInIter
                         << ", remaining Mixed " << resolve_candidates.size()
                         << " in "
                         << std::chrono::duration_cast<chrono::milliseconds>(dend-dstart).count()
                         << " ms\n";
                }

                // Any remaining Mixed leaves (resolve depth exhausted)
                // are kept: the user chose a finite resolve depth and
                // accepts that some ambiguity may persist.
                unsigned int leftoverMixed = 0;
                for (auto& q : resolve_candidates) {
                    Quadrants.push_back(std::move(q));
                    ++leftoverMixed;
                }
                // Append AllInside cells produced during resolve.
                for (auto& q : kept_non_mixed) {
                    Quadrants.push_back(std::move(q));
                }

                auto resolve_end = chrono::high_resolution_clock::now();
                cout << "    * TUSQH resolve pass done in "
                     << std::chrono::duration_cast<chrono::milliseconds>(
                            resolve_end-resolve_start).count()
                     << " ms (resolved " << totalResolved
                     << ", dropped " << totalDroppedInResolve
                     << " AllOutside, "
                     << leftoverMixed << " Mixed leaves remain)\n";
            } else {
                // No Mixed cells, nothing to do.
                Quadrants = std::move(kept_non_mixed);
            }
        }

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * windingSubdivide (TUSQH) in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count()
             << " ms" << endl;

#if (VTKOUT==true)
        {
            std::shared_ptr<FEMesh> tusqh_octree = make_shared<FEMesh>();
            saveOutputMesh(tusqh_octree, points, Quadrants, debugging);
            string tmp_name = name + "_tusqh";
            Services::WriteVTK(tmp_name, tusqh_octree);
            VolumeFractionVTKWriter::writeWindingState(tmp_name, Quadrants, points);
        }
#endif
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::computeVolumeFractions(Polyline &input, unsigned int sampleSize)
    {
        auto start_time = chrono::high_resolution_clock::now();

        WindingNumberVisitor wnv(sampleSize);
        wnv.setPolyline(&input);
        wnv.setPoints(&points);

        for (auto& q : Quadrants) {
            q.accept(&wnv);
        }

#if (VTKOUT==true)
        {
            string tmp_name = "volume_fraction_debug";
            VolumeFractionVTKWriter::writeQuadTreeWithVF(tmp_name, Quadrants, points, input);
            VolumeFractionVTKWriter::writeVFHeatmap(tmp_name, Quadrants, points);
        }
#endif

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * computeVolumeFractions (s=" << sampleSize << ") in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms" << endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::detectInsideNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        for (unsigned int i=0; i<points.size(); i++) {
            if (points[i].wasOutsideChecked()) {
                continue;
            }
            
            list<unsigned int> p_eles = points[i].getElements(), p_edges;
            points[i].outsideChecked();
            if (p_eles.empty()) {
                continue;
            }
            list<unsigned int>::const_iterator iter;
            for (iter=p_eles.begin(); iter!=p_eles.end(); ++iter) {
                const list<unsigned int> &qedges= Quadrants[*iter].getIntersectedEdges();
                //                list<unsigned int>::const_iterator qe_iter;
                //                if (qedges.empty()) {
                //                    continue;
                //                }
                // append qedges to p_edges
                p_edges.insert(p_edges.end(),qedges.begin(),qedges.end());
            }
            
            p_edges.sort();
            p_edges.unique();
            
            // p_edges: edges intersected
            if (p_edges.empty() || input.pointIsInMesh(points[i].getPoint(),p_edges)) {
                points[i].setInside();
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * detectInsideNodes in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::removeOnSurface(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        list<Quadrant> newele,removed;
        RemoveSubElementsVisitor rsv;
        rsv.setPoints(points);
        //remove elements without an inside node.
        for (unsigned int i=0; i<Quadrants.size(); i++) {
            if (Quadrants[i].isInside()) {
                newele.push_back(Quadrants[i]);
                continue;
            }
            
            if (Quadrants[i].hasIntersectedFeatures()) {
                if (Quadrants[i].accept(&rsv)) {
                    //Quadrant with feature to be removed
                    //only if average node is outside the
                    //the edges it intersects.
                    Point3D avg;
                    for (auto quaNoIdx:Quadrants[i].getPointIndex()) {
                        avg+=points[quaNoIdx].getPoint();
                    }
                    avg/=Quadrants[i].getPointIndex().size();
                    if (input.pointIsInMesh(avg,Quadrants[i].getIntersectedEdges())) {
                        newele.push_back(Quadrants[i]);
                    }
                    else {
                        removed.push_back(Quadrants[i]);
                    }
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }
            else { //FJA add a "else" here as some quadrants are inserted twice
                
                //if (Quadrants[i].removeOutsideSubElements(points)) {
                if (Quadrants[i].accept(&rsv)) {
                    removed.push_back(Quadrants[i]);
                }
                else {
                    newele.push_back(Quadrants[i]);
                }
            }
        }
        
        if (removed.empty()) {
            auto end_time = chrono::high_resolution_clock::now();
            cout << "    * RemoveOnSurface in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
            cout << " ms"<< endl;
            return;
        }
        
        //clear removed elements
        removed.clear();
        //now element std::list from Surface mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        Quadrants.assign(make_move_iterator(newele.begin()),make_move_iterator(newele.end()));
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurface in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    void Mesher::removeOnSurfaceSafe(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        list<Quadrant> newele,removed;
        
        //Keep all inside quads and all of th
        for (auto q:Quadrants) {
            if (q.isInside()) {
                newele.push_back(q);
                continue;
            }
            
            bool onein = false, feaProj = false;
            //Check inside nodes, including internal ones
            for (auto quaNoIdx:q.getSubPointIndex()) {
                if (points[quaNoIdx].isInside()) {
                    onein = true;
                    break;
                }
                if (points[quaNoIdx].isFeature()) {
                    feaProj = true;
                }
            }
            
            //if it has at least one node inside
            //we don't remove it.
            if (onein) {
                newele.push_back(q);
                continue;
            }
            
            //Compute average node of the Quadrant
            //only with corner nodes
            Point3D avg;
            for (auto quaNoIdx:q.getPointIndex()) {
                avg+=points[quaNoIdx].getPoint();
            }
            avg/=q.getPointIndex().size();
            //Sometimes the avg node is just over domain's boundary,
            //despite the fact that a (big) portion of it's still
            //inside the domain and no other Quad will cover this
            //section. For this reason the entire element is shrink
            //to 99% of it. If one node of this shrink element is still
            //inside the domain, the Octant remains. Recall that a
            //node projected onto the surface is treated as an outside
            //node.
            bool accepted = false;
            for (auto quaNoIdx:q.getPointIndex()) {
                //Point3D test = (avg + points[quaNoIdx].getPoint())/2;
                Point3D test = points[quaNoIdx].getPoint() - avg;
                test*=0.99;
                test+=avg;
                if (input.pointIsInMesh(test,q.getIntersectedEdges())) {
                    newele.push_back(q);
                    accepted = true;
                    break;
                }
            }
            if (!accepted) {
                removed.push_back(q);
            }
            
            //the most expensive case: regarding
            //the original intersected edges, does
            //this quad still have a node inside?
            //if yes we keep it. It means that all
            //of its nodes are outside or projected
            //onto the boundary so if we erase it we
            //will loose boundary representation.
            /*bool border = false;
             for (auto qp:q.getPointIndex()) {
             if (input.pointIsInMesh(points[qp].getPoint(),q.getIntersectedEdges())) {
             newele.push_back(q);
             border = true;
             break;
             }
             }
             if (!border) {
             removed.push_back(q);
             }*/
        }
        
        if (removed.empty()) {
            auto end_time = chrono::high_resolution_clock::now();
            cout << "    * RemoveOnSurfaceSafe in "
            << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
            cout << " ms"<< endl;
            return;
        }
        
        //clear removed elements
        removed.clear();
        //now element std::list from Vomule mesh can be cleared, as all remaining
        //elements are still in use and attached to newele std::list.
        Quadrants.clear();
        Quadrants.assign(std::make_move_iterator(newele.begin()),std::make_move_iterator(newele.end()));
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * RemoveOnSurfaceSafe in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::applySurfacePatterns(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //apply patterns to avoid flat, invalid and
        //poor quality elements.
        SurfaceTemplatesVisitor stv;
        stv.setPoints(points);
        stv.setPolyline(input);
        
        
        for (auto &q:Quadrants) {
            
            if (q.getPointIndex().size()!=4) {
                continue;
            }
            
            if (q.isSurface()) {
                if (!q.accept(&stv)) {
                    cout << "Error in Mesher::applySurfacePatterns: coultd't apply";
                    cout << " a surface pattern\n";
                    cout << q << "\n";
                    continue;
                }
            }
            
        }
        
        
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ApplySurfacePatterns in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    //shrink elements intersecting the envelope defined by all
    //input surfaces
    
    void Mesher::shrinkToBoundary(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the
        //surface. If after all is done, if an element counts only with "on
        //surface" and "outside" nodes, remove it.
        list<unsigned int> out_nodes;
        list<Quadrant>::iterator oiter;
        
        
        //Manage Quadrants with Features first
        for (const auto &q:Quadrants) {
            if (q.isInside()) {
                continue;
            }
            
            //Save all outside nodes and manage them after
            //dealing with all features quadrants.
            //if a node was projected to a feature
            //it will be skipped during the rest of
            //node projection.
            for (auto pIdx:q.getSubPointIndex()) {
                if (points[pIdx].isOutside()) {
                    out_nodes.push_back(pIdx);
                }
            }
            
            /*
             All features were already managed.
             
             if (!q.hasFeature()) {
             continue;
             }
             
             list<unsigned int> fs = input.getFeatureProjection(q,points);
             
             if (fs.empty()) {
             cerr << "Error at Mesher::shrinkToBoundary";
             cerr << " Quadrant labeled with feature has none\n";
             continue;
             }
             
             const vector<unsigned int> &epts = q.getPointIndex();
             
             unsigned int fsNum = fs.size(), outNo = 0;
             for (auto pIdx:epts) {
             if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
             outNo++;
             }
             }
             
             list<unsigned int>::const_iterator iter;
             
             for (iter=fs.begin(); iter!=fs.end(); ++iter) {
             
             if (outNo==0) {
             break;
             }
             
             double best = std::numeric_limits<double>::infinity();
             Point3D projected = input.getPoints()[*iter];
             unsigned int pos = 0;
             bool push = false;
             
             for (auto pIdx:epts) {
             if (points[pIdx].isOutside() && !points[pIdx].wasProjected()) {
             
             const Point3D &current = points[pIdx].getPoint();
             double dis = (current - projected).Norm();
             
             if(best>dis){
             best = dis;
             pos = pIdx;
             push = true;
             }
             }
             }
             
             if (push) {
             points[pos].setProjected();
             points[pos].setPoint(projected);
             //Feature projected flag will be used later to
             //apply surface patterns.
             points[pos].featureProjected();
             for (auto pe:points[pos].getElements()) {
             
             //this should be studied further.
             if (Quadrants.at(pe).intersectsSurface()) {
             Quadrants[pe].setSurface();
             }
             }
             outNo--;
             }
             }*/
        }
        
        //Manage non Feature Quadrants.
        out_nodes.sort();
        out_nodes.unique();
        
        for (auto p:out_nodes) {
            
            if (points.at(p).wasProjected()) {
                continue;
            }
            
            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;
            
            for (auto pe:points[p].getElements()) { //elements containing p
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }
            
            p_qInterEdges.sort();
            p_qInterEdges.unique();
            
            if (p_qInterEdges.empty()) {
                cout << "\nWarning at Mesher::shrinkToBoundary";
                cout << " no faces to project an outside node\n";
                cout << p << " n_els " << points.at(p).getElements().size() << ":";
                for (auto pe:points.at(p).getElements()) {
                    cout << " " << pe;
                }
                cout << "\n";
                continue;
            }
            
            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);
            
            //if ( current.distance(projected)< points[p].getMaxDistance() )
            {
                points[p].setPoint(projected);
                points[p].setProjected();
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ShrinkToBoundary in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;
    }
    
    
    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    
    void Mesher::projectCloseToBoundaryNodes(Polyline &input){
        auto start_time = chrono::high_resolution_clock::now();
        
        //Slow element removed (but works): from elements intersecting the
        //input domain, detect inner nodes. Project this nodes onto the surface.
        list<unsigned int> in_nodes;
        
        for (auto &q:Quadrants) {
            
            //if (!q.isSurface()) {
            if (q.isInside()) {
                continue;
            }
            //    continue;
            //}
            
            //Save the intern nodes to a list to manage them later
            //and continue the process just for feature Quadrants first.
            const vector<unsigned int> subpointindex(q.getSubPointIndex());
            for (auto pIdx:subpointindex) {
                if (points[pIdx].isInside()) {
                    in_nodes.push_back(pIdx);
                }
            }
            if (!q.hasIntersectedFeatures()) {
                continue;
            }
            
            //Manage Quadrants with Features
            list<unsigned int> fs = input.getFeatureProjection(q,points);
            
            unsigned int fsNum = fs.size();
            //we use a list and interator to erase the projected node
            //from the list of features.
            list<unsigned int>::const_iterator itFs;
            for (itFs=fs.begin(); itFs!=fs.end(); ++itFs) {
                
                //                if (fsNum==0) {
                //                    break;
                //                }
                if ( input.featureHasProjectedPt(*itFs)) {
                    break;
                }
                
                const Point3D &featProjected = input.getPoints()[*itFs];
                
                double best = std::numeric_limits<double>::infinity();
                unsigned int candidate;
                
                bool push = false;
                for (auto pIdx:subpointindex) {
                    
                    if (points[pIdx].wasProjected()) {
                        /* no more needed now we store a map Feature-ProjectedPt !!!
                         // special case: when a Feature is right on an edge/node
                         // and as already been treated in the neighboor quadrant
                         //FJA: work only if ONLY ONE feature per quad
                         // otherwise, should store the Feature number in the MeshPoint
                         // but memory consuming...
                         if (points[pIdx].isFeature()) {
                         push=false;
                         break;
                         } */
                        continue;
                    }
                    
                    const Point3D &current = points[pIdx].getPoint();
                    double dis = (current - featProjected).Norm();
                    
                    //if(points[pIdx].getMaxDistance()>dis && best>dis){
                    if(best>dis) {
                        best = dis;
                        candidate = pIdx;
                        push = true;
                    }
                }
                
                if (push) {
                    points[candidate].setProjected();
                    points[candidate].setPoint(featProjected);
                    //Feature projected flag will be used later to
                    //apply surface patterns.
                    points[candidate].setFeature();
                    input.setFeatureIndexProjectedPt(*itFs,candidate);
                    
                    for (auto pe:points[candidate].getElements()) {
                        //all the quads associated to a projected node
                        //will be labeled as surface, even if they don't
                        //intersect a boundary. This is to enable surfacePatterns
                        //for each quadrant that must be treated (note that
                        //a surface quad is not the same as an inside quad).
                        Quadrants[pe].setSurface();
                        
                        //Also, if quadrant was completely inside, then we must
                        //set as intersected edges the edges of its neighbor.
                        //This information will be used to detect if sub-elements
                        //are still inside the domain and for that, the subset
                        //of intersected edges is necessary.
                        if (Quadrants[pe].getIntersectedEdges().empty()) {
                            Quadrants[pe].setIntersectedEdges(q.getIntersectedEdges());
                        }
                    }
                    
                    fsNum--;
                }
            }
        }
        
        in_nodes.sort();
        in_nodes.unique(); //need to be sorted to call this function
        
        //move (when possible) all inner points to surface
        for (auto p:in_nodes) {
            
            if (points[p].wasProjected()) {
                //projected due to a feature => already managed.
                continue;
            }
            
            //get the faces of Quadrants sharing this node
            list<unsigned int> p_qInterEdges;
            
            //elements containing p
            for (auto pe:points.at(p).getElements()) {
                //append this to the list of edges
                p_qInterEdges.insert(p_qInterEdges.end(),
                                     Quadrants[pe].getIntersectedEdges().begin(),
                                     Quadrants[pe].getIntersectedEdges().end());
            }
            
            p_qInterEdges.sort();
            p_qInterEdges.unique();
            
            const Point3D &current = points[p].getPoint();
            Point3D projected = input.getProjection(current,p_qInterEdges);
            double dis = (current - projected).Norm();
            
            if(dis<points[p].getMaxDistance()){
                //this node has been moved to boundary, thus every element
                //sharing this node must be set as a border element in order
                //to avoid topological problems.
                points[p].setProjected();
                points[p].setPoint(projected);
                for (auto pe:points[p].getElements()) {
                    
                    //intersectsSurface() returns true when the Quad previously
                    //had non empty list of intersected edges.
                    if (Quadrants[pe].intersectsSurface()) {
                        //setSurface() will change the flags of "inside" wich
                        Quadrants[pe].setSurface();
                    }
                }
            }
        }
        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * ProjectCloseToBoundary in "
        << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count();
        cout << " ms"<< endl;

    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::computeSubcellVolumeFractions(Polyline &input,
                                               unsigned int sampleSize,
                                               double joinThreshold,
                                               const std::string &name)
    {
        auto start_time = chrono::high_resolution_clock::now();

        SubgridSampler sampler;
        unsigned int s = SubgridSampler::sanitizeSampleSize(sampleSize);

        mEdgeSubcellVF.clear();

        // Build the q_id -> vector index lookup ONCE before the edge loop.
        // EdgeInfo stores q_ids (set at Quadrant construction time) in
        // info[1]/info[2], but the Quadrants vector may be ordered
        // differently from the q_id sequence (especially after
        // resolveArchipelagos' compact-and-rebuild). We must translate
        // q_id -> index before indexing `Quadrants[]`. See BUGS_FOUND.md
        // Issue #9.
        unordered_map<unsigned int, unsigned int> qIdToIdx;
        qIdToIdx.reserve(Quadrants.size() * 2);
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            qIdToIdx[Quadrants[qi].getIndex()] = qi;
        }

        // ----- 0-cells (vertices) -----
        // For each vertex of the quadtree, build the list of incident
        // edge lengths by scanning MapEdges. Then sample s x s points
        // in the fictitious square cell centered on the vertex.
        for (unsigned int vIdx = 0; vIdx < points.size(); ++vIdx) {
            vector<pair<unsigned int, double>> incident =
                SubgridSampler::buildIncidentEdgeList(vIdx, MapEdges, points);

            SubgridSampler::Result r =
                sampler.sampleVertex(points[vIdx].getPoint(), incident, input, s);

            points[vIdx].setSubcellSampleSize(r.sampleSize);
            points[vIdx].computeSubcellVolumeFraction(r.windingNumbers);
            points[vIdx].setSubcellIsInterior(
                points[vIdx].getSubcellVolumeFraction() >= joinThreshold);
        }

        // ----- 1-cells (edges) -----
        // For each edge, get the perpendicular thickness of the
        // adjacent quads and sample s x s points in the rectangle.
        for (const auto& entry : MapEdges) {
            const QuadEdge& edge = entry.first;

            vector<double> perp =
                SubgridSampler::buildQuadPerpThickness(edge, MapEdges,
                                                       qIdToIdx, Quadrants, points);

            SubgridSampler::Result r = sampler.sampleEdge(
                points[edge[0]].getPoint(),
                points[edge[1]].getPoint(),
                perp, input, s);

            EdgeSubcellVFData data;
            data.volumeFraction = r.volumeFraction;
            data.windingNumbers = r.windingNumbers;
            data.sampleSize     = r.sampleSize;
            data.isInterior     = (r.volumeFraction >= joinThreshold);
            mEdgeSubcellVF[edge] = data;
        }

#if (VTKOUT==true)
        {
            string tmp_name = name;
            VolumeFractionVTKWriter::writeSubcellVertexVF(tmp_name, Quadrants, points,
                                                          joinThreshold);
            VolumeFractionVTKWriter::writeSubcellEdgeVF(tmp_name, Quadrants, points,
                                                        mEdgeSubcellVF, joinThreshold);
        }
#endif

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * computeSubcellVolumeFractions (s=" << s
             << ", joinThreshold=" << joinThreshold << ") in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count()
             << " ms" << endl;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    // Helper: TUSQH 1-to-5 bridge split.
    //
    // Given a quad `q` and an edge index `bridgeEdgeIdx` (0..3),
    // replaces `q` with 4 interior sub-quads (a standard 1-to-4
    // green refinement) PLUS 1 new "bridge" quad that extends
    // OUTWARD from the bridge edge by perpendicular distance
    // H/sampleSize, where H is the distance from the bridge edge to
    // the opposite edge of `q`.
    //
    // Side effects:
    //   - adds new MeshPoints (the 4 mid-edges + center from the
    //     1-to-4 split, plus 2 exterior corners for the bridge) to
    //     `points`,
    //   - adds/updates MapEdges entries for the 4 sub-quads and the
    //     bridge quad.
    //
    // Returns: vector of 5 new Quadrants (4 sub-quads + 1 bridge
    // quad) with q_ids starting at `nextQIdx`. `nextQIdx` is
    // incremented by 5 on success.
    //
    // The CALLER is responsible for:
    //   - removing `q` from `Quadrants`,
    //   - inserting the returned 5 quads in `q`'s place,
    //   - re-running computeSubcellVolumeFractions so the new edges
    //     have their fictitious cells evaluated.
    //
    // On failure (degenerate quad, SplitVisitor refuses, etc.)
    // returns an empty vector and does NOT modify any data
    // structures beyond what SplitVisitor did. The caller should
    // then leave `q` in Quadrants.
    vector<Quadrant> Mesher::bridgeSplitAtEdge(Quadrant &q,
                                                Quadrant &q2,
                                                unsigned int bridgeEdgeIdx,
                                                unsigned int sampleSize,
                                                unsigned int &nextQIdx,
                                                bool doManifoldSplit)
    {
        (void)doManifoldSplit;  // Option B-subdiv is planned but not
                                 // implemented yet (kept in the API
                                 // for future use). Option A is the
                                 // current behaviour.

        vector<Quadrant> result;

        if (bridgeEdgeIdx > 3) return result;

        const vector<unsigned int> &pi = q.getPointIndex();
        if (pi.size() != 4) return result;

        // q2 is the neighbour on the other side of the bridge edge.
        // MapEdges guarantees info[2] is a valid quad; we don't
        // re-validate its corners because the cubical complex may
        // already have refined or split it but kept its id slot for
        // bookkeeping (in that case the BFS still considers it a
        // neighbour via the MapEdges entry). The caller has already
        // confirmed both sides are interior and in different
        // components via the candidate filter.
        (void)q2;

        if (sampleSize == 0) return result;

        // Identify the bridge edge endpoints and the opposite edge
        // endpoints. The bridge edge is edge i = (pi[i], pi[(i+1)%4]).
        // The opposite edge is (pi[(i+2)%4], pi[(i+3)%4]).
        const unsigned int e0 = pi[bridgeEdgeIdx];
        const unsigned int e1 = pi[(bridgeEdgeIdx + 1) % 4];
        const unsigned int e2 = pi[(bridgeEdgeIdx + 2) % 4];
        const unsigned int e3 = pi[(bridgeEdgeIdx + 3) % 4];

        const Point3D &p_e0 = points[e0].getPoint();
        const Point3D &p_e1 = points[e1].getPoint();
        const Point3D &p_e2 = points[e2].getPoint();
        const Point3D &p_e3 = points[e3].getPoint();

        // Issue #2 (TUSQH §3.4): outward direction via quad centroid.
        //
        // The original implementation computed the outward direction as
        // `-normalize(mid_opposite - mid_bridge)`, which assumes a
        // convex, axis-aligned quad. For rotated or non-convex quads
        // the midpoint of the opposite edge does not necessarily lie
        // on the "interior" side. Using the centroid (mean of the 4
        // corners) is robust: it always lies inside the convex hull.
        //
        // H is still defined as the centroid-to-edge distance, so the
        // bridge thickness H/sampleSize is consistent with the paper's
        // 1-to-5 template geometry.
        Point3D mid_bridge = (p_e0 + p_e1) * 0.5;
        Point3D dir_exterior_unit = computeExteriorDirection(
            q, bridgeEdgeIdx, points);
        if (dir_exterior_unit.Norm() < 1e-12) return result;  // degenerate

        Point3D centroid = (p_e0 + p_e1 + p_e2 + p_e3) * 0.25;
        double H = (centroid - mid_bridge).Norm();
        if (H < 1e-12) return result;

        double bridge_thickness = H / (double)sampleSize;

        Point3D p_e0_ext = p_e0 + dir_exterior_unit * bridge_thickness;
        Point3D p_e1_ext = p_e1 + dir_exterior_unit * bridge_thickness;

        // Apply the existing 1-to-4 SplitVisitor. It writes:
        //   - 4 sub-quads into `new_eles` (each a vector of 4 point
        //     indices that already reference future indices in
        //     `points`)
        //   - 4 mid-edge points + 1 center point into `new_pts`
        //     (in that order; some may be missing if their edge was
        //     already split)
        //   - new sub-edges and updated midpoint info into MapEdges.
        SplitVisitor sv;
        list<Point3D> new_pts;
        vector<vector<unsigned int>> new_eles;
        vector<vector<Point3D>> clipping_coords;
        vector<Quadrant> processed;
        map<unsigned int, unsigned int> idx_pos_map;
        list<pair<unsigned int, unsigned int>> to_balance;

        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setNewEles(new_eles);
        sv.setClipping(clipping_coords);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(to_balance);
        sv.setStartIndex(nextQIdx);

        bool visit_ok = sv.visit(&q);
        if (!visit_ok || new_eles.size() != 4) {
            return result;
        }

        // Append SplitVisitor's new points (mid-edges and center)
        // to `points`. After this, the point indices inside
        // `new_eles[i]` are valid in the global `points` vector.
        for (const auto& p : new_pts) {
            points.emplace_back(MeshPoint(p));
        }

        // Add the 2 new exterior corners for the bridge quad.
        unsigned int e0_ext_idx = (unsigned int)points.size();
        points.emplace_back(MeshPoint(p_e0_ext));
        unsigned int e1_ext_idx = (unsigned int)points.size();
        points.emplace_back(MeshPoint(p_e1_ext));

        // Find the midpoint on the bridge edge (which SplitVisitor
        // stored as `info[0]` of the (e0, e1) edge entry).
        unsigned int m_bridge;
        {
            unsigned int a = e0, b = e1;
            if (a > b) std::swap(a, b);
            QuadEdge ke(a, b);
            auto it = MapEdges.find(ke);
            if (it == MapEdges.end() || it->second[0] == 0) {
                // SplitVisitor didn't split the bridge edge (it was
                // already split by an earlier operation). For now
                // abort; we don't support nested splits.
                return result;
            }
            m_bridge = it->second[0];
        }

        // q_id for the bridge quad: 4 sub-quads get nextQIdx..+3,
        // bridge gets nextQIdx+4.
        unsigned int bridge_q_id = nextQIdx + 4;

        // Update the 2 half-edges of the bridge edge so that the
        // bridge quad appears as the second quad (info[2]).
        auto updateHalfEdge = [&](unsigned int p1, unsigned int p2) {
            unsigned int a = p1, b = p2;
            if (a > b) std::swap(a, b);
            QuadEdge ke(a, b);
            auto it = MapEdges.find(ke);
            if (it != MapEdges.end() && it->second[2] == std::numeric_limits<unsigned int>::max()) {
                it->second[2] = bridge_q_id;
            }
        };
        updateHalfEdge(e0, m_bridge);
        updateHalfEdge(m_bridge, e1);

        // Add 3 new edges for the bridge quad (the side opposite
        // the shared bridge edge).
        auto addBridgeEdge = [&](unsigned int p1, unsigned int p2) {
            unsigned int a = p1, b = p2;
            if (a > b) std::swap(a, b);
            QuadEdge ke(a, b, true);
            auto it = MapEdges.find(ke);
            if (it == MapEdges.end()) {
                MapEdges.emplace(ke, EdgeInfo(0, bridge_q_id,
                    std::numeric_limits<unsigned int>::max()));
            } else if (it->second[2] == std::numeric_limits<unsigned int>::max()) {
                it->second[2] = bridge_q_id;
            } else if (it->second[1] == std::numeric_limits<unsigned int>::max()) {
                it->second[1] = bridge_q_id;
            }
        };
        addBridgeEdge(e1, e1_ext_idx);
        addBridgeEdge(e1_ext_idx, e0_ext_idx);
        addBridgeEdge(e0_ext_idx, e0);

        // Build the 5 new Quadrants. The 4 sub-quads come from
        // `new_eles`; the bridge quad has corners (e0, e1, e1_ext,
        // e0_ext) traversed CCW.
        unsigned short ref_level = (unsigned short)(q.getRefinementLevel() + 1);

        for (unsigned int i = 0; i < 4; ++i) {
            Quadrant sq(new_eles[i], ref_level, nextQIdx++);
            result.push_back(std::move(sq));
        }

        vector<unsigned int> bridge_corners = {e0, e1, e1_ext_idx, e0_ext_idx};
        Quadrant bq(bridge_corners, ref_level, nextQIdx++);
        result.push_back(std::move(bq));

        // Issue #8 (TUSQH §3.4, manifoldness): instead of erasing the
        // (e0, e1) entry from MapEdges (which made the bridge quad
        // topologically isolated), UPDATE it so that:
        //   - info[0] (midpoint) stays at m_bridge (set by SplitVisitor)
        //   - info[1] becomes the bridge quad's q_id (was stale: pointed
        //     to the removed quad q)
        //   - info[2] becomes the neighbour q2's q_id (was already set
        //     by SplitVisitor but we re-assert it for safety)
        //
        // This way the bridge quad is reachable from q2 in BFS via
        // the full (e0, e1) edge, so the bridge quad joins q2's
        // component. The bridge quad geometrically lies INSIDE q2's
        // rectangular fictitious cell (the bridge edge is the shared
        // boundary), so this matches TUSQH §3.4's "bridge quad is on
        // the side of the neighbour with smaller sum-of-VFs".
        //
        // Note: this is Option A (topology fix only). The bridge quad
        // still geometrically overlaps q2's rectangular cell, so
        // strict manifoldness is not achieved. Option B-subdiv
        // (splitting q2 and discarding 2 of its sub-quads) would fix
        // manifoldness but is not yet implemented.
        {
            unsigned int a = e0, b = e1;
            if (a > b) std::swap(a, b);
            QuadEdge ke(a, b);
            auto it = MapEdges.find(ke);
            if (it == MapEdges.end()) {
                // Should not happen: SplitVisitor created this entry.
                MapEdges.emplace(ke, EdgeInfo(m_bridge, bridge_q_id,
                    q2.getIndex()));
            } else {
                it->second[0] = m_bridge;
                it->second[1] = bridge_q_id;
                it->second[2] = q2.getIndex();
            }
        }

        return result;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    Point3D Mesher::computeExteriorDirection(const Quadrant &q,
                                              unsigned int bridgeEdgeIdx,
                                              const vector<MeshPoint> &points) const
    {
        Point3D zero(0.0, 0.0, 0.0);
        const vector<unsigned int> &pi = q.getPointIndex();
        if (bridgeEdgeIdx > 3 || pi.size() != 4) return zero;

        const unsigned int e0 = pi[bridgeEdgeIdx];
        const unsigned int e1 = pi[(bridgeEdgeIdx + 1) % 4];

        const Point3D &p_e0 = points[e0].getPoint();
        const Point3D &p_e1 = points[e1].getPoint();
        Point3D mid_bridge = (p_e0 + p_e1) * 0.5;

        // Quad centroid = mean of the 4 corners. The interior
        // direction is from the bridge edge midpoint to the centroid;
        // the exterior direction is its negation, normalised.
        Point3D centroid(0.0, 0.0, 0.0);
        for (unsigned int i = 0; i < 4; ++i) {
            centroid += points[pi[i]].getPoint();
        }
        centroid *= 0.25;

        Point3D dir_interior = centroid - mid_bridge;
        double H = dir_interior.Norm();
        if (H < 1e-12) return zero;  // degenerate quad

        return (-dir_interior) / H;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    // Paper-faithful: an "interior" cell of the cubical complex is one
    // whose volume fraction is positive (AllInside) or whose winding
    // numbers are ambiguous (Mixed, sits on the boundary). AllOutside
    // cells are kept in Quadrants so MapEdges retains the topology
    // between components, but they are skipped everywhere a mesh
    // entity is needed.
    bool Mesher::isInteriorCell(const Quadrant &q) {
        switch (q.getWindingState()) {
            case WindingState::AllInside:
            case WindingState::Mixed:
                return true;
            case WindingState::AllOutside:
            case WindingState::Unknown:
            default:
                return false;
        }
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    bool Mesher::isEdgeOnDomainBoundary(const QuadEdge &edge,
                                        const Quadrant &adjacentQuad,
                                        const Polyline &input,
                                        unsigned int sampleSize,
                                        double joinThreshold) const
    {
        // We need to know which edge index inside `adjacentQuad`
        // corresponds to `edge`. The adjacentQuad is the quad that
        // owns the side `info[1]` of this edge; `edge` is one of its
        // four edges.
        const vector<unsigned int> &pi = adjacentQuad.getPointIndex();
        unsigned int edgeIdx = 4;
        for (unsigned int i = 0; i < pi.size(); ++i) {
            unsigned int a = pi[i];
            unsigned int b = pi[(i + 1) % pi.size()];
            if ((a == edge[0] && b == edge[1]) ||
                (a == edge[1] && b == edge[0])) {
                edgeIdx = i;
                break;
            }
        }
        if (edgeIdx == 4) return false;  // not actually adjacent

        unsigned int s = SubgridSampler::sanitizeSampleSize(sampleSize);
        if (s == 0) return false;

        // Exterior direction (outward from the quad through this edge).
        Point3D dir_ext_unit = computeExteriorDirection(
            adjacentQuad, edgeIdx, points);
        if (dir_ext_unit.Norm() < 1e-12) return false;

        // Edge axis (along the edge).
        const Point3D &p_a = points[edge[0]].getPoint();
        const Point3D &p_b = points[edge[1]].getPoint();
        double dx = p_b[0] - p_a[0];
        double dy = p_b[1] - p_a[1];
        double edgeLen = std::sqrt(dx * dx + dy * dy);
        if (edgeLen <= 0.0) return false;
        double ax = dx / edgeLen;
        double ay = dy / edgeLen;

        // Perpendicular thickness: use the quad perpendicular height H.
        Point3D centroid(0.0, 0.0, 0.0);
        for (unsigned int i = 0; i < pi.size(); ++i) {
            centroid += points[pi[i]].getPoint();
        }
        centroid *= 0.25;
        Point3D mid_edge = (p_a + p_b) * 0.5;
        double H = (centroid - mid_edge).Norm();
        if (H <= 0.0) H = edgeLen;
        if (H <= 0.0) return false;

        const double halfAlong = 0.5 * edgeLen;
        const double stepAlong = edgeLen / static_cast<double>(s);
        const double stepPerp  = H     / static_cast<double>(s);

        // Sample s x s points in the fictitious rectangle of the edge,
        // restricted to the +dir_ext half-plane (j index goes from 0
        // to s-1, all positive offset along dir_ext).
        double sum = 0.0;
        unsigned int n = 0;
        for (unsigned int i = 0; i < s; ++i) {
            for (unsigned int j = 0; j < s; ++j) {
                double u = -halfAlong + (i + 0.5) * stepAlong;
                double v = (j + 0.5) * stepPerp;  // >= 0
                double ex = dir_ext_unit[0];
                double ey = dir_ext_unit[1];
                double x = mid_edge[0] + u * ax + v * ex;
                double y = mid_edge[1] + u * ay + v * ey;
                Point3D sample(x, y, 0.0);
                int wn = input.windingNumber(sample);
                sum += static_cast<double>(wn);
                ++n;
            }
        }
        if (n == 0) return false;

        double vf_exterior = sum / static_cast<double>(n);

        // Strict semantics: vf_exterior < joinThreshold ⇔ the edge's
        // exterior is mostly outside the domain ⇔ the edge sits on
        // the global domain boundary.
        return vf_exterior < joinThreshold;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    void Mesher::resolveArchipelagos(Polyline &input,
                                      unsigned int sampleSize,
                                      double joinThreshold,
                                      unsigned int minComponentCells,
                                      const string& name,bool Aliasing)
    {
        auto start_time = chrono::high_resolution_clock::now();

        // ---- 0) Bridge-joining loop ----
        // TUSQH paper §3.4 step 3: "For any pair of connected
        // components, if the edges that connect them are interior
        // to the geometry, the components are joined using templates
        // along those edges."
        //
        // Paper-faithful implementation (feature/paper-faithful-bridge):
        //
        // The cubical complex is the full background grid kept in
        // `Quadrants` (see Step 1 of this branch: AllOutside cells
        // are no longer dropped). AllOutside cells are invisible to
        // the BFS (`isInteriorCell` returns false) but they keep the
        // edges in MapEdges intact, so the algorithm can find pairs
        // of interior components whose separation is a chain of
        // exterior cells.
        //
        // Each iteration:
        //   A) BFS over interior cells (AllInside + Mixed) using the
        //      preserved MapEdges topology. This assigns a component
        //      label to every interior quad.
        //   B) For each edge in MapEdges where BOTH sides are interior
        //      quads in DIFFERENT components, with interior sub-cell
        //      VF (>= joinThreshold) and NOT on the global domain
        //      boundary (Issue #1), perform a 1-to-5 split on the
        //      quad on side info[1]. The bridge quad extends OUTWARD
        //      from info[1]'s centroid toward info[2]'s centroid,
        //      which is on the other side of the shared edge — the
        //      bridge quad lands ON info[2]'s territory, connecting
        //      the two components.
        //   C) Re-compute sub-cell VFs so the new bridge quads'
        //      outer edges can be evaluated in the next iteration.
        //
        // Repeat until no new bridges are added or maxBridgeIterations
        // is reached. The iteration cap prevents runaway growth on
        // pathological inputs.
        const int maxBridgeIterations = 5;
        int totalBridgesAdded = 0;

        // q_id counter for new bridge quads. Start at one past the
        // largest existing q_id so we don't collide.
        unsigned int nextQIdx = 0;
        for (auto& q : Quadrants) {
            unsigned int qi = q.getIndex();
            if (qi != std::numeric_limits<unsigned int>::max() && qi + 1 > nextQIdx) {
                nextQIdx = qi + 1;
            }
        }

        for (int bridgeIter = 0; bridgeIter < maxBridgeIterations; ++bridgeIter) {
            // Refresh sub-cell VFs so the new edges from the previous
            // iteration can be considered. (No-op on iter 0 since
            // the caller has already computed them.)
            if (bridgeIter > 0) {
                computeSubcellVolumeFractions(input, sampleSize, joinThreshold);
            }

            // Build q_id -> vector index lookup for the current
            // Quadrants vector.
            std::unordered_map<unsigned int, unsigned int> qIdToIdx;
            qIdToIdx.reserve(Quadrants.size() * 2);
            for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
                qIdToIdx[Quadrants[qi].getIndex()] = qi;
            }

            // ---- A) BFS over the FULL cubical complex ----
            // Per the paper §3.4, components are detected over the
            // background grid (which includes AllOutside cells).
            // Treating AllOutside cells as BFS nodes lets the BFS hop
            // across chains of exterior cells and identify true
            // connected regions. The BFS labels EVERY cell, so two
            // interior cells separated only by an exterior corridor
            // end up in the SAME component — which is the correct
            // paper-faithful interpretation.
            vector<int> compOfQuad(Quadrants.size(), -1);
            int numComponents = 0;
            for (unsigned int start = 0; start < Quadrants.size(); ++start) {
                if (compOfQuad[start] != -1) continue;
                std::queue<unsigned int> q;
                q.push(start);
                compOfQuad[start] = numComponents;
                while (!q.empty()) {
                    unsigned int qi = q.front(); q.pop();
                    const auto& pi = Quadrants[qi].getPointIndex();
                    for (unsigned int e = 0; e < pi.size(); ++e) {
                        unsigned int a = pi[e];
                        unsigned int b = pi[(e + 1) % pi.size()];
                        if (a > b) std::swap(a, b);
                        QuadEdge ke(a, b);  // default ctor sorts
                        auto it = MapEdges.find(ke);
                        if (it == MapEdges.end()) continue;
                        const EdgeInfo& info = it->second;
                        for (unsigned int k = 1; k <= 2; ++k) {
                            unsigned int otherQId = info[k];
                            if (otherQId == std::numeric_limits<unsigned int>::max()) continue;
                            auto itMap = qIdToIdx.find(otherQId);
                            if (itMap == qIdToIdx.end()) continue;
                            unsigned int otherIdx = itMap->second;
                            if (compOfQuad[otherIdx] != -1) continue;
                            compOfQuad[otherIdx] = numComponents;
                            q.push(otherIdx);
                        }
                    }
                }
                ++numComponents;
            }

            // ---- B) Bridge candidate collection ----
            // For each MapEdge whose BOTH sides are interior quads in
            // DIFFERENT components, with interior sub-cell VF and NOT
            // on the global domain boundary, add a 1-to-5 split. The
            // bridge quad extends OUTWARD from info[1]'s centroid,
            // across the shared edge into info[2]'s territory — this
            // is what makes the two components join.
            vector<QuadEdge> bridgeEdges;
            int filteredBoundary = 0;

            // Nuevo --------------------------------
            bool requiresBalancing = false;
            vector<unsigned int> quadsToRefine;
            // ----------------------------------------

            for (const auto& entry : MapEdges) {
                const EdgeInfo& info = entry.second;
                if (info[1] == std::numeric_limits<unsigned int>::max()) continue;
                if (info[2] == std::numeric_limits<unsigned int>::max()) continue;
                auto it1 = qIdToIdx.find(info[1]);
                auto it2 = qIdToIdx.find(info[2]);
                if (it1 == qIdToIdx.end() || it2 == qIdToIdx.end()) continue;
                unsigned int qi1 = it1->second;
                unsigned int qi2 = it2->second;

                if (!isInteriorCell(Quadrants[qi1])) continue;
                if (!isInteriorCell(Quadrants[qi2])) continue;
                if (compOfQuad[qi1] == compOfQuad[qi2]) continue;

                // Nuevo ---------------------------------------------
                int level1 = (int)Quadrants[qi1].getRefinementLevel();
                int level2 = (int)Quadrants[qi2].getRefinementLevel();

                if (level1 != level2) {
                    requiresBalancing = true;
                    // Guardamos la celda más grande (menor nivel) para subdividirla
                    if (level1 < level2) quadsToRefine.push_back(qi1);
                    else quadsToRefine.push_back(qi2);

                    // Bloqueamos la creación del puente temporalmente.
                    // Continuamos evaluando otras aristas para juntar todas
                    // las celdas que necesiten subdivisión en esta pasada.
                    continue;
                }
                // ----------------------------------------------------------

                // The edge's fictitious cell must be interior (>= threshold).
                auto vfIt = mEdgeSubcellVF.find(entry.first);
                if (vfIt == mEdgeSubcellVF.end()) continue;
                if (vfIt->second.volumeFraction < joinThreshold) continue;
                // Issue #1: edge must NOT be on the global domain
                // boundary (the other side is genuinely outside the
                // geometry). In that case adding a bridge would land
                // a quad in true exterior.
                if (isEdgeOnDomainBoundary(entry.first, Quadrants[qi1],
                                           input, sampleSize, joinThreshold)) {
                    ++filteredBoundary;
                    continue;
                }
                bridgeEdges.push_back(entry.first);
            }

            // Nuevo ----------------------------------------------------------

            if (requiresBalancing) {
                // Deduplicar: varias aristas pueden pedir refinar la misma celda.
                std::sort(quadsToRefine.begin(), quadsToRefine.end());
                quadsToRefine.erase(std::unique(quadsToRefine.begin(), quadsToRefine.end()), quadsToRefine.end());

                cout << "    [bridge iter " << bridgeIter
                     << "] Equilibrado JIT: subdividiendo " << quadsToRefine.size()
                     << " celdas dispares (1-irregularidad forzada).\n";

                // Ordenar en orden DECRECIENTE de índice para que el push_back
                // al final del vector (dentro de balanceSplitQuad) no invalide
                // los índices aún por procesar en este lote.
                std::sort(quadsToRefine.begin(), quadsToRefine.end(),
                          std::greater<unsigned int>());

                unsigned int balancedCount = 0;
                unsigned int totalSurvivors = 0;
                for (unsigned int qi : quadsToRefine) {
                    if (qi >= Quadrants.size()) continue;
                    auto newOnes = balanceSplitQuad(input, qi, sampleSize, nextQIdx);
                    if (!newOnes.empty()) {
                        ++balancedCount;
                        totalSurvivors += (unsigned int)newOnes.size();
                    }
                }

                cout << "    [bridge iter " << bridgeIter
                     << "] Equilibrado JIT completado: " << balancedCount
                     << "/" << quadsToRefine.size()
                     << " subdivisiones aplicadas ("
                     << totalSurvivors << " sub-quads retenidos).\n";

                // La iteración se gastó en equilibrado; el BFS del próximo
                // ciclo se ejecutará sobre la malla ya balanceada. Saltamos
                // la creación de puentes de esta pasada. El equilibrado es
                // "gratis": reservamos las maxBridgeIterations iteraciones
                // reales para construcción de puentes, evitando agotar el
                // presupuesto en subdivisiones JIT consecutivas.
                --bridgeIter;
                continue;
            }

            // ---------------------------------------------------------------

            if (filteredBoundary > 0) {
                cout << "    [bridge iter " << bridgeIter
                     << "] filtered " << filteredBoundary
                     << " boundary edges (domain boundary)\n";
            }

            if (bridgeEdges.empty()) break;

            // AQUI SE INSERTAN LOS TEMPLATES
            int bridgesThisIter = 0;
            for (const auto& ke : bridgeEdges) {
                auto meIt = MapEdges.find(ke);
                if (meIt == MapEdges.end()) continue;
                unsigned int q_id = meIt->second[1];
                unsigned int q2_id = meIt->second[2];
                if (q_id == std::numeric_limits<unsigned int>::max()) continue;
                if (q2_id == std::numeric_limits<unsigned int>::max()) continue;

                auto qiIt = qIdToIdx.find(q_id);
                auto q2iIt = qIdToIdx.find(q2_id);
                if (qiIt == qIdToIdx.end() || q2iIt == qIdToIdx.end()) continue;
                unsigned int qi = qiIt->second;
                unsigned int q2i = q2iIt->second;

                // Locate the edge index in the quad's corner list.
                const auto& pi = Quadrants[qi].getPointIndex();
                unsigned int bridgeEdgeIdx = 4;
                for (unsigned int i = 0; i < 4; ++i) {
                    unsigned int p1 = pi[i];
                    unsigned int p2 = pi[(i + 1) % 4];
                    if ((p1 == ke[0] && p2 == ke[1]) ||
                        (p1 == ke[1] && p2 == ke[0])) {
                        bridgeEdgeIdx = i;
                        break;
                    }
                }
                if (bridgeEdgeIdx == 4) continue;

                // Perform the 1-to-5 split. Pass q2 so the bridge quad
                // can be topologically registered as a neighbour of q2
                // (Issue #8 fix: MapEdges update instead of erase).
                vector<Quadrant> newQuads = bridgeSplitAtEdge(
                    Quadrants[qi], Quadrants[q2i], bridgeEdgeIdx,
                    sampleSize, nextQIdx, false);
                if (newQuads.empty()) continue;

                // Replace the original quad with the first new quad
                // and append the rest.
                Quadrants[qi] = std::move(newQuads[0]);
                for (unsigned int k = 1; k < newQuads.size(); ++k) {
                    Quadrants.push_back(std::move(newQuads[k]));
                }

                // Classify every new quad so its WindingState is set
                // (the Quadrant constructor leaves it as Unknown).
                // This includes the first sub-quad at position qi, which
                // replaced the original quad and would otherwise stay
                // Unknown forever — that would make isInteriorCell()
                // return false in the final drop filter and discard an
                // otherwise-interior cell. The original code only
                // classified the 4 appended quads and missed the
                // replacement whenever qi was not the last position.
                //
                // WindingNumberVisitor::visit only computes the per-sample
                // winding numbers + VF; it does NOT set WindingState.
                // WindingState must be derived from the winding numbers
                // we just computed: all wn>0 → AllInside, all wn==0 →
                // AllOutside, otherwise → Mixed. (Equivalent to what
                // WindingNumberSubdivisionVisitor does, but without the
                // subdivision side-effect.)
                auto classifyFromWindingNumbers = [](Quadrant &q) {
                    const auto &wns = q.getWindingNumbers();
                    bool anyPos = false, anyZero = false;
                    for (double wn : wns) {
                        if (wn > 0.0) anyPos = true;
                        else anyZero = true;
                    }
                    if (anyPos && !anyZero)
                        q.setWindingState(WindingState::AllInside);
                    else if (!anyPos && anyZero)
                        q.setWindingState(WindingState::AllOutside);
                    else
                        q.setWindingState(WindingState::Mixed);
                };
                if (!Quadrants.empty()) {
                    WindingNumberVisitor wnv(sampleSize);
                    wnv.setPolyline(&input);
                    wnv.setPoints(&points);
                    Quadrants[qi].accept(&wnv);
                    classifyFromWindingNumbers(Quadrants[qi]);
                    const unsigned int startIdx =
                        Quadrants.size() - (unsigned int)newQuads.size() + 1;
                    for (unsigned int k = startIdx; k < Quadrants.size(); ++k) {
                        Quadrants[k].accept(&wnv);
                        classifyFromWindingNumbers(Quadrants[k]);
                    }
                }
                ++bridgesThisIter;
                // TEMPLATES
            }

            if (bridgesThisIter == 0) break;
            totalBridgesAdded += bridgesThisIter;

            cout << "    [bridge iter " << bridgeIter
                 << "] added " << bridgesThisIter
                 << " bridges (total=" << totalBridgesAdded << ")\n";

#if (VTKOUT==true)
            {
                // Step 5' — Per-iteration bridge state. Snapshot of the
                // quadtree right after the bridges of this iteration
                // have been added (BEFORE the next iteration's BFS
                // relabels components). Useful to see how the bridge
                // quads physically extend outward and which components
                // are merging.
                std::shared_ptr<FEMesh> iter_mesh = make_shared<FEMesh>();
                if (!Quadrants.empty()) {
                    saveOutputMesh(iter_mesh, points, Quadrants);
                }
                string tmp_name = name + "_bridge_iter" + std::to_string(bridgeIter);
                Services::WriteVTK(tmp_name, iter_mesh);
            }
#endif
        }

#if (VTKOUT==true)
        {
            // Step 5.5a — Post-bridges winding-state snapshot. Same
            // geometry as output_postbridges.vtk but with per-cell
            // WindingState (Unknown=0 / AllInside=1 / AllOutside=2 /
            // Mixed=3) as CELL_DATA, plus volume_fraction and
            // refinement_level for context.
            //
            // Diagnostic target: confirm the Fix A in resolveArchipelagos
            // (Mesher.cpp:3649-3668) classified the first sub-quad after
            // bridgeSplitAtEdge; if any cell shows Unknown=0 here, the
            // BFS will mark it unvisited and the drop filter will discard
            // it. Also useful to inspect bridge quads that may have been
            // placed inside what should be AllInside regions (those would
            // appear as Mixed=3 cells surrounded by AllInside=1, which
            // is a strong hint of an issue in bridgeSplitAtEdge or
            // computeExteriorDirection for non-AABB quads).
            string tmp_name = name + "_postbridge";
            if (!Quadrants.empty()) {
                VolumeFractionVTKWriter::writeWindingState(
                    tmp_name, Quadrants, points);
                unsigned int nUnknown = 0, nAllInside = 0,
                             nAllOutside = 0, nMixed = 0;
                for (const auto& q : Quadrants) {
                    switch (q.getWindingState()) {
                        case WindingState::Unknown:    ++nUnknown; break;
                        case WindingState::AllInside:  ++nAllInside; break;
                        case WindingState::AllOutside: ++nAllOutside; break;
                        case WindingState::Mixed:      ++nMixed; break;
                    }
                }
                cout << "    [postbridge winding] Unknown=" << nUnknown
                     << " AllInside=" << nAllInside
                     << " AllOutside=" << nAllOutside
                     << " Mixed=" << nMixed << "\n";
            }
        }
#endif

#if (VTKOUT==true)
        {
            // Step 5.5 — Post-bridges snapshot. State of the quadtree
            // after the bridge loop has finished but BEFORE the final
            // BFS relabels components and BEFORE the small-component
            // drop. Diff against output_bridge_iter<N>.vtk shows only
            // the final-iter bridge additions; diff against
            // output_postbfs.vtk shows the BFS relabel effect; diff
            // against output_predrop.vtk shows what the drop removes.
            string tmp_name = name + "_postbridges";
            if (!Quadrants.empty()) {
                std::shared_ptr<FEMesh> mesh = make_shared<FEMesh>();
                saveOutputMesh(mesh, points, Quadrants);
                Services::WriteVTK(tmp_name, mesh);
                cout << "  Wrote: " << tmp_name << ".vtk ("
                     << Quadrants.size() << " cells, "
                     << totalBridgesAdded << " bridges added)\n";
            } else {
                cout << "  Skipped: " << tmp_name
                     << ".vtk (Quadrants empty)\n";
            }
        }
#endif

        // ---- 1) Find connected components of the quadtree ----
        // Final pass (used for the small-component drop below and for
        // the "numComponents" reported in the summary). Same logic as
        // the inner BFS above: full cubical complex BFS, all cells
        // are labeled.
        vector<int> compOfQuad(Quadrants.size(), -1);
        std::unordered_map<unsigned int, unsigned int> qIdToIdx;
        qIdToIdx.reserve(Quadrants.size() * 2);
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            qIdToIdx[Quadrants[qi].getIndex()] = qi;
        }
        int numComponents = 0;

        for (unsigned int start = 0; start < Quadrants.size(); ++start) {
            if (compOfQuad[start] != -1) continue;
            // BFS over quads via shared edges, full cubical complex
            // (both interior and AllOutside cells).
            std::queue<unsigned int> q;
            q.push(start);
            compOfQuad[start] = numComponents;
            while (!q.empty()) {
                unsigned int qi = q.front(); q.pop();
                const auto& pi = Quadrants[qi].getPointIndex();
                for (unsigned int e = 0; e < pi.size(); ++e) {
                    unsigned int a = pi[e];
                    unsigned int b = pi[(e + 1) % pi.size()];
                    // MapEdges stores the canonical (sorted) key, so
                    // we must sort a/b before constructing the key.
                    if (a > b) std::swap(a, b);
                    QuadEdge ke(a, b);  // default ctor sorts
                    auto it = MapEdges.find(ke);
                    if (it == MapEdges.end()) continue;
                    const EdgeInfo& info = it->second;
                    for (unsigned int k = 1; k <= 2; ++k) {
                        unsigned int otherQId = info[k];
                        if (otherQId == std::numeric_limits<unsigned int>::max()) continue;
                        auto itMap = qIdToIdx.find(otherQId);
                        if (itMap == qIdToIdx.end()) continue;
                        unsigned int otherIdx = itMap->second;
                        if (compOfQuad[otherIdx] != -1) continue;
                        compOfQuad[otherIdx] = numComponents;
                        q.push(otherIdx);
                    }
                }
            }
            ++numComponents;
        }

#if (VTKOUT==true)
        {
            // Step 5.6 — Post-BFS snapshot. Same geometry as
            // output_postbridges.vtk, but each cell now carries its
            // component label as a CELL_DATA scalar. Colour by
            // `component_id` in ParaView to see which cells belong to
            // the same connected component.
            string tmp_name = name + "_postbfs";
            if (!Quadrants.empty()) {
                vector<double> compVals(Quadrants.size(), -1.0);
                for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
                    if (compOfQuad[qi] >= 0) {
                        compVals[qi] = static_cast<double>(compOfQuad[qi]);
                    }
                }
                VolumeFractionVTKWriter::writeQuadTreeWithCellArray(
                    tmp_name, Quadrants, points, "component_id", compVals);
                cout << "  Wrote: " << tmp_name << ".vtk ("
                     << numComponents << " components)\n";
            } else {
                cout << "  Skipped: " << tmp_name
                     << ".vtk (Quadrants empty)\n";
            }
        }
#endif

        // ---- 2) Drop "small" components ----
        // Per the paper §3.4: "The remaining connected components that
        // contain fewer than a user-defined number of highest-
        // dimensional cells are removed." Components are computed
        // over the full cubical complex (interior + AllOutside),
        // so the size is the total number of cells in the component.
        // Cells that are NOT interior (AllOutside) are dropped from
        // the output regardless of component membership, since they
        // were only preserved for the BFS / bridge topology.
        vector<int> compSize(numComponents, 0);
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            if (compOfQuad[qi] >= 0) compSize[compOfQuad[qi]]++;
        }

        int droppedQuads = 0;
vector<bool> keepQuad(Quadrants.size(), false);
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            // Keep: interior cell, in a visited component, component
            // size >= min. AllOutside cells in kept components are
            // dropped from the output regardless of component membership, since they
            // were preserved only for the BFS / bridge topology.
            if (compOfQuad[qi] >= 0 &&
                compSize[compOfQuad[qi]] >= (int)minComponentCells &&
                isInteriorCell(Quadrants[qi])) {
                keepQuad[qi] = true;
            }
        }

#if (VTKOUT==true)
        {
            // Step 5.7 — Pre-drop snapshot. Same geometry as
            // output_postbfs.vtk but each cell now has a 0/1 flag
            // `will_keep` indicating whether it survives the small-
            // component drop. Colour by `will_keep` in ParaView to
            // see exactly which cells are about to be removed.
            string tmp_name = name + "_predrop";
            if (!Quadrants.empty()) {
                vector<double> keepVals(Quadrants.size(), 0.0);
                unsigned int willKeepCount = 0;
                for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
                    if (keepQuad[qi]) {
                        keepVals[qi] = 1.0;
                        ++willKeepCount;
                    }
                }
                VolumeFractionVTKWriter::writeQuadTreeWithCellArray(
                    tmp_name, Quadrants, points, "will_keep", keepVals);
                cout << "  Wrote: " << tmp_name << ".vtk ("
                     << willKeepCount << " will keep, "
                     << Quadrants.size() - willKeepCount << " will drop)\n";
            } else {
                cout << "  Skipped: " << tmp_name
                     << ".vtk (Quadrants empty)\n";
            }
        }
#endif

        // Compact the Quadrants vector and rebuild MapEdges from kept
        // quads. We deliberately re-use EdgeVisitor so the resulting
        // edge map matches what the rest of the pipeline expects.
        vector<Quadrant> kept;
        kept.reserve(Quadrants.size());
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            if (keepQuad[qi]) kept.push_back(Quadrants[qi]);
        }
        droppedQuads = (int)(Quadrants.size() - kept.size());

        // Rebuild MapEdges from scratch using the kept quads. We have
        // to renumber q_ids so that EdgeInfo can still reference them
        // by index. The simplest stable mapping is "old index -> new
        // index" if the quad was kept.
        for (unsigned int oldIdx = 0; oldIdx < Quadrants.size(); ++oldIdx) {
            if (keepQuad[oldIdx]) {
                unsigned int newIdx = (unsigned int)kept.size();
                // Quadrant::q_id is set at construction time. We
                // bypass immutability by relying on the fact that the
                // pipeline only uses q_id for VTK/debug, while
                // adjacency is by vector position. So we just leave
                // q_id as-is.
                (void)newIdx;
            }
        }

        // Rebuild MapEdges.
        map<QuadEdge, EdgeInfo> newMapEdges;
        for (unsigned int newIdx = 0; newIdx < kept.size(); ++newIdx) {
            const auto& pi = kept[newIdx].getPointIndex();
            for (unsigned int e = 0; e < pi.size(); ++e) {
                unsigned int a = pi[e];
                unsigned int b = pi[(e + 1) % pi.size()];
                QuadEdge ke(a, b, true);
                auto it = newMapEdges.find(ke);
                if (it == newMapEdges.end()) {
                    newMapEdges.emplace(ke, EdgeInfo(0, newIdx,
                        std::numeric_limits<unsigned int>::max()));
                } else {
                    it->second[2] = newIdx;
                }
            }
        }
        MapEdges = std::move(newMapEdges);
        Quadrants = std::move(kept);

        // Rebuild points[i].elements for the (now smaller) Quadrants
        // vector. linkElementsToNodes() will be called by the
        // pipeline downstream.
        for (auto& v : points) v.clearElements();
        for (unsigned int qi = 0; qi < Quadrants.size(); ++qi) {
            for (unsigned int vIdx : Quadrants[qi].getPointIndex()) {
                if (vIdx < points.size()) points[vIdx].addElement(qi);
            }
        }

        // ---- 3) Drop edge sub-cell VF entries whose quads vanished ----
        vector<QuadEdge> stale;
        for (const auto& entry : mEdgeSubcellVF) {
            const QuadEdge& ke = entry.first;
            auto it = MapEdges.find(ke);
            if (it == MapEdges.end()) stale.push_back(ke);
        }
        for (const auto& ke : stale) mEdgeSubcellVF.erase(ke);

#if (VTKOUT==true)
        {
            string tmp_name = name + "_postarchipelago";
            std::shared_ptr<FEMesh> post_mesh = make_shared<FEMesh>();
            // Guard: saveOutputMesh accesses tmp_Quadrants[0] (line
            // 1623) without an empty check, so it would segfault if
            // all components were dropped by -L. This is a pre-existing
            // issue in the mesher that becomes reachable now that the
            // Issue #1 filter can leave Quadrants empty.
            if (!Quadrants.empty()) {
                saveOutputMesh(post_mesh, points, Quadrants);
                Services::WriteVTK(tmp_name, post_mesh);
                cout << "  Wrote: " << tmp_name << ".vtk\n";
            } else {
                cout << "  Skipped: " << tmp_name
                     << ".vtk (all components dropped, Quadrants empty)\n";
            }
        }
#endif

        auto end_time = chrono::high_resolution_clock::now();
        cout << "    * resolveArchipelagos: " << numComponents
             << " components, " << totalBridgesAdded << " bridges added, "
             << droppedQuads
             << " small-component quads dropped (min=" << minComponentCells
             << ") in "
             << std::chrono::duration_cast<chrono::milliseconds>(end_time-start_time).count()
             << " ms" << endl;

        //------------ Desarollo Joaquin Fierro ---------------------
        if (Aliasing == true){  
            list<Quadrant> tmp_Quadrants ;
            tmp_Quadrants.assign(make_move_iterator(Quadrants.begin()),make_move_iterator(Quadrants.end()));
            unsigned int new_q_idx = tmp_Quadrants.size();
            cout << "Se ha activado el uso de templates \n"<< endl;
            QuadAliasing qa;
            qa.setQuadrant(tmp_Quadrants);
            qa.setPoints(points);
            qa.setActualIndex(new_q_idx);
            cout<< "Detectando Pinches" << endl;
            qa.getPinches();
            cout<< "Generando Pinches" << endl;
            qa.CreateTemplates();

            std::shared_ptr<FEMesh> templates_octree=make_shared<FEMesh>();
            saveOutputMesh(templates_octree,points,tmp_Quadrants);
            string tmp_name = name + "_templatesPinches";
            Services::WriteVTK(tmp_name,templates_octree);
        } 
        //------------ Fin desarrollo -------------------------

        
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    // JIT balancing helper: 1-to-4 split on Quadrants[qi] with full
    // bookkeeping (MapEdges rebuild, AllOutside filter, WindingState
    // classification). Mirrors the bridge loop's pattern (Mesher.cpp
    // ~3721-3776) and splitQuadrants' pattern (~625-704).
    //
    // On success, replaces Quadrants[qi] with the first surviving
    // sub-quad and pushes the remaining survivors at the end of
    // Quadrants. Updates MapEdges for the new sub-quads (SplitVisitor
    // splits the outer edges into halves and creates the internal
    // cross edges). Appends the new MeshPoints created by SplitVisitor
    // (4 mid-edges + 1 center) to `points`. Bumps nextQIdx past the
    // assigned q_ids.
    //
    // AllOutside children are discarded: their MapEdges entries get
    // info[1]/info[2] set to numeric_limits::max via
    // EdgeVisitor::removeEdges, matching splitQuadrants' filter.
    //
    // Returns the survivors (1 to 4). Returns an empty vector if the
    // quad should not be split (out of range, visitor failed, or all
    // children ended up AllOutside).
    //
    // IMPORTANT: the caller MUST process the list of indices to split
    // in DECREASING order so that push_back at the end of Quadrants
    // does not invalidate subsequent indices in the batch.
    //
    // MapEdges coherence: SplitVisitor mutates MapEdges in place
    // (splitting the 4 outer edges at their midpoints and adding 4
    // cross-edges to the center). For AllOutside children we call
    // EdgeVisitor::removeEdges to set their side pointers to max so
    // the MapEdges entries are not orphaned. If ALL children end up
    // AllOutside, the function reverts every MapEdges change made by
    // SplitVisitor (restores info[0]=0 on the 4 original edges and
    // erases the 8 half-edges + 4 cross-edges) so the BFS of the
    // next bridge iteration cannot walk phantom topology. This is
    // the same protocol used by bridgeSplitAtEdge (see lines 3109+).
    std::vector<Quadrant> Mesher::balanceSplitQuad(
        Polyline &input,
        unsigned int qi,
        unsigned int sampleSize,
        unsigned int &nextQIdx)
    {
        std::vector<Quadrant> result;
        if (qi >= Quadrants.size()) return result;

        Quadrant quad = Quadrants[qi];

        unsigned short qrl = quad.getRefinementLevel();
        list<unsigned int> &inter_edges = quad.getIntersectedEdges();

        // Reuse SplitVisitor just like splitQuadrants and bridgeSplitAtEdge.
        SplitVisitor sv;
        list<Point3D> new_pts;
        vector<vector<unsigned int>> split_elements;
        vector<vector<Point3D>> clipping_coords;
        vector<Quadrant> processed;
        map<unsigned int, unsigned int> idx_pos_map;
        list<pair<unsigned int, unsigned int>> to_balance;

        sv.setPoints(points);
        sv.setMapEdges(MapEdges);
        sv.setNewPts(new_pts);
        sv.setNewEles(split_elements);
        sv.setClipping(clipping_coords);
        sv.setProcessedQuadVector(processed);
        sv.setMapProcessed(idx_pos_map);
        sv.setToBalanceList(to_balance);
        sv.setStartIndex(nextQIdx);

        if (!quad.accept(&sv) || split_elements.size() != 4) {
            return result;
        }

        // Append SplitVisitor's new MeshPoints (mid-edges + center).
        if (!new_pts.empty()) {
            points.reserve(points.size() + new_pts.size());
            points.insert(points.end(), new_pts.begin(), new_pts.end());
        }

        // Build the 4 candidate children. We track survivors in
        // `survivors` and discarded indices in `discardIndices` so we
        // can later clean up MapEdges for the discarded ones.
        std::vector<Quadrant> survivors;
        std::vector<unsigned int> discardIndices;
        const unsigned int qidBase = nextQIdx;
        for (unsigned int j = 0; j < split_elements.size(); ++j) {
            unsigned int childQId = qidBase + j;
            Quadrant o(split_elements[j],
                       static_cast<unsigned short>(qrl + 1),
                       childQId);

            bool keep = false;
            if (inter_edges.empty()) {
                keep = true;
            } else {
                IntersectionsVisitor iv(true);
                iv.setPolyline(input);
                iv.setEdges(inter_edges);
                iv.setCoords(clipping_coords[j]);
                if (o.accept(&iv)) {
                    keep = true;
                } else if (isItIn(input, inter_edges, clipping_coords[j])) {
                    keep = true;
                }
            }

            if (keep) {
                survivors.push_back(std::move(o));
            } else {
                // AllOutside: clean up its MapEdges entries (set the
                // side pointers to max so the entry is not orphaned).
                EdgeVisitor::removeEdges(&o, MapEdges);
                discardIndices.push_back(j);
            }
        }

        if (survivors.empty()) {
            // All children ended up AllOutside. SplitVisitor already
            // mutated MapEdges: it stamped info[0] on the 4 original
            // outer edges with their midpoints, and inserted 8 new
            // half-edges plus 4 cross-edges pointing at q_ids in
            // [qidBase, qidBase+3]. Since we are NOT inserting any
            // Quadrant with those q_ids, we must revert MapEdges so
            // the BFS of the next iteration does not walk phantom
            // topology. The discarded-children side pointers have
            // already been zeroed by EdgeVisitor::removeEdges above.
            const vector<unsigned int> &pi = quad.getPointIndex();
            unsigned int mid01 = MapEdges[QuadEdge(pi[0], pi[1])][0];
            unsigned int mid12 = MapEdges[QuadEdge(pi[1], pi[2])][0];
            unsigned int mid23 = MapEdges[QuadEdge(pi[2], pi[3])][0];
            unsigned int mid30 = MapEdges[QuadEdge(pi[3], pi[0])][0];
            // The center is one corner common to all 4 sub-quads; for
            // the first child it is the 3rd corner (split_elements[0][2]).
            const unsigned int center = split_elements[0][2];

            MapEdges[QuadEdge(pi[0], pi[1])][0] = 0;
            MapEdges[QuadEdge(pi[1], pi[2])][0] = 0;
            MapEdges[QuadEdge(pi[2], pi[3])][0] = 0;
            MapEdges[QuadEdge(pi[3], pi[0])][0] = 0;

            MapEdges.erase(QuadEdge(pi[0],  mid01));
            MapEdges.erase(QuadEdge(mid01,   pi[1]));
            MapEdges.erase(QuadEdge(pi[1],  mid12));
            MapEdges.erase(QuadEdge(mid12,   pi[2]));
            MapEdges.erase(QuadEdge(pi[2],  mid23));
            MapEdges.erase(QuadEdge(mid23,   pi[3]));
            MapEdges.erase(QuadEdge(pi[3],  mid30));
            MapEdges.erase(QuadEdge(mid30,   pi[0]));

            MapEdges.erase(QuadEdge(mid01,  center));
            MapEdges.erase(QuadEdge(mid12,  center));
            MapEdges.erase(QuadEdge(mid23,  center));
            MapEdges.erase(QuadEdge(mid30,  center));

            nextQIdx = qidBase;
            return result;
        }

        // Bump nextQIdx past all 4 children regardless of survival
        // (q_ids were assigned at construction; keeping them ensures
        // no q_id reuse, even for discarded AllOutside children).
        nextQIdx = qidBase + 4;

        // Replace Quadrants[qi] with the first survivor and push the
        // rest at the end of the vector.
        Quadrants[qi] = std::move(survivors[0]);
        for (unsigned int k = 1; k < survivors.size(); ++k) {
            Quadrants.push_back(std::move(survivors[k]));
        }

        // Classify WindingState for each survivor (Unknown is the
        // Quadrant constructor default; without classification,
        // isInteriorCell would return false and the new sub-quads
        // would be discarded by the final drop filter).
        auto classifyFromWindingNumbers = [](Quadrant &q) {
            const auto &wns = q.getWindingNumbers();
            bool anyPos = false, anyZero = false;
            for (double wn : wns) {
                if (wn > 0.0) anyPos = true;
                else anyZero = true;
            }
            if (anyPos && !anyZero)
                q.setWindingState(WindingState::AllInside);
            else if (!anyPos && anyZero)
                q.setWindingState(WindingState::AllOutside);
            else
                q.setWindingState(WindingState::Mixed);
        };

        WindingNumberVisitor wnv(sampleSize);
        wnv.setPolyline(&input);
        wnv.setPoints(&points);
        // Re-visit Quadrants[qi] (which now holds the first survivor)
        // and the appended survivors at the end of the vector.
        Quadrants[qi].accept(&wnv);
        classifyFromWindingNumbers(Quadrants[qi]);
        const unsigned int startIdx =
            Quadrants.size() - survivors.size() + 1;
        for (unsigned int k = startIdx; k < Quadrants.size(); ++k) {
            Quadrants[k].accept(&wnv);
            classifyFromWindingNumbers(Quadrants[k]);
        }

        // Collect the survivors as currently stored in Quadrants for
        // the caller (so the caller has stable copies).
        std::vector<Quadrant> returned;
        returned.reserve(survivors.size());
        returned.push_back(Quadrants[qi]);
        for (unsigned int k = startIdx; k < Quadrants.size(); ++k) {
            returned.push_back(Quadrants[k]);
        }
        return returned;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------

    // Heurística para el error de cuadrícula gruesa de TUSQH. Si el quadtree comienza
    // con `coarseThreshold` o menos celdas raíz y la polilínea tiene
    // al menos `minSegmentsForTrigger` segmentos, el criterio de bobinado de TUSQH
    // puede submuestrear las celdas iniciales gigantes (ej. con el valor por
    // defecto -N 2 = 4 muestras por celda, una celda de 50x50 se clasifica como
    // TodoAdentro o TodoAfuera y nunca se subdivide). Se pre-refina
    // uniformemente hasta el nivel `baseLevel` usando la función heredada `generateQuadtreeMesh`
    // para que TUSQH tenga suficientes celdas pequeñas y sus cálculos sean significativos.
    //
    // El pre-refinamiento está limitado: solo se activa cuando la cuadrícula
    // inicial es muy gruesa (por defecto: <= 2 celdas), la polilínea es
    // lo suficientemente compleja (por defecto: >= 100 segmentos; esto excluye casos
    // de regresión simples como unit_square.poly donde TUSQH produce
    // legítimamente una sola celda), y respeta la profundidad máxima solicitada `maxDepth`
    // (limitada a baseLevel, para que el bucle de TUSQH todavía tenga margen para
    // subdividir más si es necesario).
    bool Mesher::preRefineForTusqh(Polyline &input,
                                   unsigned int maxDepth,
                                   list<RefinementRegion *> &all_reg,
                                   const string &name,
                                   const bool &debugging,
                                   bool Aliasing,
                                   unsigned short coarseThreshold,
                                   unsigned short baseLevel,
                                   unsigned int minSegmentsForTrigger) {
        if (Quadrants.size() > coarseThreshold) {
            return false;
        }
        if (input.getEdges().size() < minSegmentsForTrigger) {
            return false;
        }
        unsigned short effectiveBaseLevel = baseLevel;
        if (maxDepth < effectiveBaseLevel) {
            effectiveBaseLevel = static_cast<unsigned short>(maxDepth);
        }
        if (effectiveBaseLevel == 0) {
            return false;
        }

        cout << "    * preRefineForTusqh: initial grid has "
             << Quadrants.size() << " cell(s) (<= " << coarseThreshold
             << ") and input has " << input.getEdges().size()
             << " segments (>= " << minSegmentsForTrigger
             << "); pre-refining uniformly to level " << effectiveBaseLevel
             << " to avoid TUSQH winding undersampling\n";

        unsigned int pre_q_idx = Quadrants.size();
        generateQuadtreeMesh(effectiveBaseLevel, input, all_reg, name,
                             0, effectiveBaseLevel, debugging,
                             pre_q_idx, Aliasing);

        cout << "    * preRefineForTusqh: quadtree now has "
             << Quadrants.size() << " cells before TUSQH subdivision\n";

        return true;
    }
}
