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
 * @file VolumeFractionVTKWriter.h
 * @author Claudio Lobos, Fabrice Jaillet
 * @version 0.1
 * @brief Debug output for volume fraction visualization
 **/

#ifndef VolumeFractionVTKWriter_h
#define VolumeFractionVTKWriter_h

#include <string>
#include <vector>

namespace Clobscode
{
    class Quadrant;
    class MeshPoint;
    class Polyline;

    class VolumeFractionVTKWriter
    {
    public:

        static bool writeQuadTreeWithVF(const std::string& name,
                                        const std::vector<Quadrant>& quadrants,
                                        const std::vector<MeshPoint>& points,
                                        const Polyline& geo);

        static bool writeVFHeatmap(const std::string& name,
                                   const std::vector<Quadrant>& quadrants,
                                   const std::vector<MeshPoint>& points);

        // Debug visualisation of the TUSQH winding-state classification:
        //   0 = Unknown
        //   1 = AllInside
        //   2 = AllOutside
        //   3 = Mixed
        static bool writeWindingState(const std::string& name,
                                      const std::vector<Quadrant>& quadrants,
                                      const std::vector<MeshPoint>& points);
    };
}

#endif