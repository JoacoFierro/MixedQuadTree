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
  along with this program.  If not, see <http://www.gnu.org/licenses/lgpl.txt>.
*/
/**
* @file RunStats.h
* @brief POD struct that aggregates metrics from a single mesher run.
*        Used to feed Services::AppendStatsCSV.
**/

#ifndef RunStats_h
#define RunStats_h 1

#include <string>
#include <vector>

namespace Clobscode
{

struct RunStats
{
    // --- Identificación ---
    std::string in_name;
    std::string out_name;
    std::string timestamp;

    // --- Quadtree / subdivisión ---
    unsigned int n_initial_quadrants = 0;
    unsigned int n_final_quadrants   = 0;

    // --- TUSQH ---
    unsigned int tusqh_iter_count        = 0;
    unsigned int preserved_outside       = 0;
    unsigned int mixed_before_resolve    = 0;
    unsigned int leftover_mixed          = 0;
    unsigned int tusqh_refined_count     = 0;

    // --- Archipiélagos ---
    unsigned int num_components           = 0;
    unsigned int total_bridges_added      = 0;
    unsigned int dropped_quads            = 0;
    unsigned int min_component_cells      = 0;
    unsigned int n_postbridge_unknown     = 0;
    unsigned int n_postbridge_allinside   = 0;
    unsigned int n_postbridge_alloutside  = 0;
    unsigned int n_postbridge_mixed       = 0;

    // --- Aliasing (pinches) ---
    unsigned int n_pinches_detected       = 0;
    unsigned int n_pinch_templates        = 0;
    unsigned int n_archipelago_templates  = 0;
    double       template_pct_quads       = 0.0;

    // --- Mesh final (derivado de FEMesh) ---
    unsigned int n_elements_total         = 0;
    unsigned int n_quads                  = 0;
    unsigned int n_triangles              = 0;
    double       quad_triangle_ratio      = 0.0;
    unsigned int n_points                 = 0;
    unsigned int n_outside_nodes          = 0;

    // --- Calidad (de FEMesh decoration) ---
    double min_angle_quad = 0.0;
    double max_angle_quad = 0.0;
    double min_angle_tri  = 0.0;
    double max_angle_tri  = 0.0;

    // --- Distribución por nivel (de FEMesh::ref_levels) ---
    unsigned int max_ref_level           = 0;
    unsigned int min_ref_level           = 0;
    double       mean_ref_level          = 0.0;
    static const unsigned int MAX_LEVEL_SLOTS = 16;
    unsigned int n_elements_at_level[MAX_LEVEL_SLOTS] = {0};

    // --- Tiempos (ms) ---
    double gen_time_ms   = 0.0;
    double write_time_ms = 0.0;
    double all_time_ms   = 0.0;

    // Returns the canonical 55-column header line (no trailing newline).
    static const std::vector<std::string> &columnNames();

    // Returns the corresponding cell values formatted as strings.
    // Strings that may contain commas / quotes / newlines are quoted.
    std::vector<std::string> toRow() const;
};

}

#endif
