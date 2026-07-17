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
  along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>.
*/
/**
  * @file SubcellVFData.h
  * @brief POD structs used to carry per-cell subgrid-sampled volume
  *        fraction data (TUSQH paper, §3.3) through Mesher and the
  *        VTK writer without modifying the on-disk EdgeInfo layout.
  **/

#ifndef SubcellVFData_h
#define SubcellVFData_h 1

#include <vector>

namespace Clobscode
{
    /// Subgrid-sampled volume fraction data for a quadtree edge.
    /// Held in a side-map keyed by QuadEdge so it does not interfere
    /// with the binary format of EdgeInfo.
    struct EdgeSubcellVFData {
        double volumeFraction;
        std::vector<double> windingNumbers;
        unsigned int sampleSize;
        bool isInterior;  // VF >= joinThreshold
    };
}

#endif
