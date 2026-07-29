/*
  <Mix-mesher: region type. This program generates a mixed-elements 2D mesh>

  Copyright (C) <2013,2018>  <Claudio Lobos> All rights reserved.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
*/

#include "RunStats.h"

#include <cstdio>
#include <cmath>
#include <sstream>
#include <iomanip>

namespace Clobscode
{

static std::string escapeCsv(const std::string &s)
{
    bool needs_quote = false;
    for (char c : s) {
        if (c == ',' || c == '"' || c == '\n' || c == '\r') {
            needs_quote = true;
            break;
        }
    }
    if (!needs_quote) return s;

    std::string out;
    out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) {
        if (c == '"') out.push_back('"');
        out.push_back(c);
    }
    out.push_back('"');
    return out;
}

static std::string fmtUnsigned(unsigned int v)
{
    char buf[32];
    std::snprintf(buf, sizeof(buf), "%u", v);
    return std::string(buf);
}

static std::string fmtDouble(double v)
{
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%.6g", v);
    return std::string(buf);
}

const std::vector<std::string> &RunStats::columnNames()
{
    static const std::vector<std::string> names = {
        // Identificación
        "timestamp", "in_name", "out_name",
        // Flags
        "ref_level", "decoration", "mSampleSize",
        "useTusqh", "tusqhSampleSize", "refineOnEdgeIntersect", "tusqhExtraResolveDepth",
        "useSubgrid", "subgridSampleSize", "subgridJoinThreshold", "subgridMinComponentCells",
        "pinchDetectionMode", "Aliasing",
        // Quadtree / TUSQH
        "n_initial_quadrants", "n_final_quadrants",
        "tusqh_iter_count", "preserved_outside", "mixed_before_resolve",
        "leftover_mixed", "tusqh_refined_count",
        // Archipiélagos
        "num_components", "total_bridges_added", "dropped_quads", "min_component_cells",
        "n_postbridge_unknown", "n_postbridge_allinside",
        "n_postbridge_alloutside", "n_postbridge_mixed",
        // Aliasing
        "n_pinches_detected", "n_pinch_templates", "n_archipelago_templates",
        "template_pct_quads",
        // Mesh
        "n_elements_total", "n_quads", "n_triangles", "quad_triangle_ratio",
        "n_points", "n_outside_nodes",
        // Calidad
        "min_angle_quad", "max_angle_quad", "min_angle_tri", "max_angle_tri",
        // Distribución por nivel
        "max_ref_level", "min_ref_level", "mean_ref_level",
        "n_elements_at_level_0",  "n_elements_at_level_1",
        "n_elements_at_level_2",  "n_elements_at_level_3",
        "n_elements_at_level_4",  "n_elements_at_level_5",
        "n_elements_at_level_6",  "n_elements_at_level_7",
        "n_elements_at_level_8",  "n_elements_at_level_9",
        "n_elements_at_level_10", "n_elements_at_level_11",
        "n_elements_at_level_12", "n_elements_at_level_13",
        "n_elements_at_level_14", "n_elements_at_level_15",
        // Tiempos
        "gen_time_ms", "write_time_ms", "all_time_ms"
    };
    return names;
}

std::vector<std::string> RunStats::toRow() const
{
    std::vector<std::string> row;
    row.reserve(columnNames().size());

    // Identificación
    row.push_back(escapeCsv(timestamp));
    row.push_back(escapeCsv(in_name));
    row.push_back(escapeCsv(out_name));

    // Flags - empty placeholders; populated by Main.cpp from the local
    // flag variables, since Mesher doesn't see them directly.  We use
    // std::numeric_limits<unsigned int>::max() as a sentinel meaning
    // "not set" for unsigneds, and an empty string for booleans.
    row.push_back("");  // ref_level
    row.push_back("");  // decoration
    row.push_back("");  // mSampleSize
    row.push_back("");  // useTusqh
    row.push_back("");  // tusqhSampleSize
    row.push_back("");  // refineOnEdgeIntersect
    row.push_back("");  // tusqhExtraResolveDepth
    row.push_back("");  // useSubgrid
    row.push_back("");  // subgridSampleSize
    row.push_back("");  // subgridJoinThreshold
    row.push_back("");  // subgridMinComponentCells
    row.push_back("");  // pinchDetectionMode
    row.push_back("");  // Aliasing

    // Quadtree / TUSQH
    row.push_back(fmtUnsigned(n_initial_quadrants));
    row.push_back(fmtUnsigned(n_final_quadrants));
    row.push_back(fmtUnsigned(tusqh_iter_count));
    row.push_back(fmtUnsigned(preserved_outside));
    row.push_back(fmtUnsigned(mixed_before_resolve));
    row.push_back(fmtUnsigned(leftover_mixed));
    row.push_back(fmtUnsigned(tusqh_refined_count));

    // Archipiélagos
    row.push_back(fmtUnsigned(num_components));
    row.push_back(fmtUnsigned(total_bridges_added));
    row.push_back(fmtUnsigned(dropped_quads));
    row.push_back(fmtUnsigned(min_component_cells));
    row.push_back(fmtUnsigned(n_postbridge_unknown));
    row.push_back(fmtUnsigned(n_postbridge_allinside));
    row.push_back(fmtUnsigned(n_postbridge_alloutside));
    row.push_back(fmtUnsigned(n_postbridge_mixed));

    // Aliasing
    row.push_back(fmtUnsigned(n_pinches_detected));
    row.push_back(fmtUnsigned(n_pinch_templates));
    row.push_back(fmtUnsigned(n_archipelago_templates));
    row.push_back(fmtDouble(template_pct_quads));

    // Mesh
    row.push_back(fmtUnsigned(n_elements_total));
    row.push_back(fmtUnsigned(n_quads));
    row.push_back(fmtUnsigned(n_triangles));
    row.push_back(fmtDouble(quad_triangle_ratio));
    row.push_back(fmtUnsigned(n_points));
    row.push_back(fmtUnsigned(n_outside_nodes));

    // Calidad
    row.push_back(fmtDouble(min_angle_quad));
    row.push_back(fmtDouble(max_angle_quad));
    row.push_back(fmtDouble(min_angle_tri));
    row.push_back(fmtDouble(max_angle_tri));

    // Distribución por nivel
    row.push_back(fmtUnsigned(max_ref_level));
    row.push_back(fmtUnsigned(min_ref_level));
    row.push_back(fmtDouble(mean_ref_level));
    for (unsigned int i = 0; i < MAX_LEVEL_SLOTS; ++i) {
        row.push_back(fmtUnsigned(n_elements_at_level[i]));
    }

    // Tiempos
    row.push_back(fmtDouble(gen_time_ms));
    row.push_back(fmtDouble(write_time_ms));
    row.push_back(fmtDouble(all_time_ms));

    return row;
}

}
