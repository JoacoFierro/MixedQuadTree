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
 * @file VolumeFractionVTKWriter.cpp
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief Debug output for volume fraction visualization
 **/

#include "VolumeFractionVTKWriter.h"
#include "../Quadrant.h"
#include "../MeshPoint.h"
#include "../Polyline.h"
#include "../Point3D.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::vector;
using std::string;

namespace Clobscode
{

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool VolumeFractionVTKWriter::writeQuadTreeWithVF(const std::string& name,
                                                      const std::vector<Quadrant>& quadrants,
                                                      const std::vector<MeshPoint>& points,
                                                      const Polyline& geo)
    {
        if (quadrants.empty()) {
            std::cerr << "No quadrants to write\n";
            return false;
        }

        string vol_name = name + ".vtk";

        FILE* f = fopen(vol_name.c_str(), "wt");
        if (!f) {
            std::cerr << "Cannot open file: " << vol_name << "\n";
            return false;
        }

        // Header VTK
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "Volume Fraction QuadTree Debug Output\n");
        fprintf(f, "ASCII\n\n");

        // Los puntos son los nodos del quadtree ( MeshPoints )
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(f, "POINTS %u float\n", (unsigned int)points.size());

        for (unsigned int i = 0; i < points.size(); i++) {
            const Point3D& p = points[i].getPoint();
            fprintf(f, "%+1.8E %+1.8E %+1.8E\n", p[0], p[1], p[2]);
        }

        // Contar celdas (cada quadrant es una celda quad)
        unsigned int num_cells = quadrants.size();
        unsigned int connectivity = 0;
        for (const auto& q : quadrants) {
            const auto& idx = q.getPointIndex();
            connectivity += idx.size() + 1;
        }

        fprintf(f, "\nCELLS %u %u\n", num_cells, connectivity);

        for (const auto& q : quadrants) {
            const auto& idx = q.getPointIndex();
            fprintf(f, "%u", (unsigned int)idx.size());
            for (auto i : idx) {
                fprintf(f, " %u", i);
            }
            fprintf(f, "\n");
        }

        // CELL_TYPES
        fprintf(f, "\nCELL_TYPES %u\n", num_cells);
        for (unsigned int i = 0; i < num_cells; i++) {
            if (i % 30 == 0) fprintf(f, "\n");
            fprintf(f, "9 "); // VTK_QUAD
        }

        // CELL_DATA con volume fractions
        fprintf(f, "\n\nCELL_DATA %u\n", num_cells);
        fprintf(f, "SCALARS volume_fraction double 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");

        for (unsigned int i = 0; i < quadrants.size(); i++) {
            if (i % 30 == 0) fprintf(f, "\n");
            fprintf(f, "%+1.8E\n", quadrants[i].getVolumeFraction());
        }

        // Agregar también información de si intersecta superficie
        fprintf(f, "\nSCALARS intersects_surface int 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (const auto& q : quadrants) {
            fprintf(f, "%d\n", q.intersectsSurface() ? 1 : 0);
        }

        // Agregar winding numbers como vector si está disponible
        fprintf(f, "\nSCALARS refinement_level unsigned_int 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (const auto& q : quadrants) {
            fprintf(f, "%u\n", q.getRefinementLevel());
        }

        fclose(f);
        std::cout << "  Wrote: " << vol_name << "\n";
        return true;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool VolumeFractionVTKWriter::writeVFHeatmap(const std::string& name,
                                                 const std::vector<Quadrant>& quadrants,
                                                 const std::vector<MeshPoint>& points)
    {
        // Esta función escribe los centroides de sub_elements con sus winding numbers
        // para visualizar el campo de volume fraction como un heatmap

        if (quadrants.empty()) {
            std::cerr << "No quadrants to write\n";
            return false;
        }

        string vol_name = name + "_samples.vtk";

        FILE* f = fopen(vol_name.c_str(), "wt");
        if (!f) {
            std::cerr << "Cannot open file: " << vol_name << "\n";
            return false;
        }

        // Header VTK
        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "s×s Grid Sample Points with Winding Numbers Debug Output\n");
        fprintf(f, "ASCII\n\n");

        // Recolectar puntos de la grilla s×s de cada cuadrante
        vector<Point3D> centroids;
        vector<double> wn_values;

        for (const auto& q : quadrants) {
            unsigned int s = q.getSampleSize();
            const auto& wn = q.getWindingNumbers();

            for (unsigned int i = 0; i < s; ++i) {
                for (unsigned int j = 0; j < s; ++j) {
                    Point3D sample = q.getSamplePoint(i, j, points);
                    centroids.push_back(sample);
                    unsigned int idx = i * s + j;
                    wn_values.push_back(wn[idx]);
                }
            }
        }

        unsigned int num_points = centroids.size();

        // Escribir puntos
        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(f, "POINTS %u float\n", num_points);

        for (const auto& p : centroids) {
            fprintf(f, "%+1.8E %+1.8E %+1.8E\n", p[0], p[1], p[2]);
        }

        // Celdas - cada punto es un vértice
        fprintf(f, "\nCELLS %u %u\n", num_points, num_points * 2);
        for (unsigned int i = 0; i < num_points; i++) {
            fprintf(f, "1 %u\n", i);
        }

        // CELL_TYPES - todos vertices
        fprintf(f, "\nCELL_TYPES %u\n", num_points);
        for (unsigned int i = 0; i < num_points; i++) {
            fprintf(f, "1 "); // VTK_VERTEX
        }

        // POINT_DATA con winding numbers
        fprintf(f, "\nPOINT_DATA %u\n", num_points);
        fprintf(f, "SCALARS winding_number double 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (double v : wn_values) {
            fprintf(f, "%+1.8E\n", v);
        }

        fclose(f);
        std::cout << "  Wrote: " << vol_name << "\n";
        return true;
    }

    //--------------------------------------------------------------------------------
    //--------------------------------------------------------------------------------
    bool VolumeFractionVTKWriter::writeWindingState(const std::string& name,
                                                    const std::vector<Quadrant>& quadrants,
                                                    const std::vector<MeshPoint>& points)
    {
        if (quadrants.empty()) {
            std::cerr << "No quadrants to write\n";
            return false;
        }

        string vol_name = name + "_winding_state.vtk";

        FILE* f = fopen(vol_name.c_str(), "wt");
        if (!f) {
            std::cerr << "Cannot open file: " << vol_name << "\n";
            return false;
        }

        fprintf(f, "# vtk DataFile Version 2.0\n");
        fprintf(f, "TUSQH Winding State Classification Debug Output\n");
        fprintf(f, "ASCII\n\n");

        fprintf(f, "DATASET UNSTRUCTURED_GRID\n");
        fprintf(f, "POINTS %u float\n", (unsigned int)points.size());
        for (unsigned int i = 0; i < points.size(); i++) {
            const Point3D& p = points[i].getPoint();
            fprintf(f, "%+1.8E %+1.8E %+1.8E\n", p[0], p[1], p[2]);
        }

        unsigned int num_cells = quadrants.size();
        unsigned int connectivity = 0;
        for (const auto& q : quadrants) {
            const auto& idx = q.getPointIndex();
            connectivity += idx.size() + 1;
        }
        fprintf(f, "\nCELLS %u %u\n", num_cells, connectivity);
        for (const auto& q : quadrants) {
            const auto& idx = q.getPointIndex();
            fprintf(f, "%u", (unsigned int)idx.size());
            for (auto i : idx) {
                fprintf(f, " %u", i);
            }
            fprintf(f, "\n");
        }

        fprintf(f, "\nCELL_TYPES %u\n", num_cells);
        for (unsigned int i = 0; i < num_cells; i++) {
            if (i % 30 == 0) fprintf(f, "\n");
            fprintf(f, "9 ");
        }

        fprintf(f, "\n\nCELL_DATA %u\n", num_cells);
        fprintf(f, "SCALARS winding_state int 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (const auto& q : quadrants) {
            int s = 0;
            switch (q.getWindingState()) {
                case WindingState::Unknown:    s = 0; break;
                case WindingState::AllInside:  s = 1; break;
                case WindingState::AllOutside: s = 2; break;
                case WindingState::Mixed:      s = 3; break;
            }
            fprintf(f, "%d\n", s);
        }

        // Also include the volume fraction for context.
        fprintf(f, "\nSCALARS volume_fraction double 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (const auto& q : quadrants) {
            fprintf(f, "%+1.8E\n", q.getVolumeFraction());
        }

        fprintf(f, "\nSCALARS refinement_level unsigned_int 1\n");
        fprintf(f, "LOOKUP_TABLE default\n");
        for (const auto& q : quadrants) {
            fprintf(f, "%u\n", q.getRefinementLevel());
        }

        fclose(f);
        std::cout << "  Wrote: " << vol_name << "\n";
        return true;
    }

}