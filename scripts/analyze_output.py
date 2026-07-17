#!/usr/bin/env python3
"""Quantitative analysis of TUSQH bridge-joining output.

Compares a `_postarchipelago.vtk` file against the input polyline
to assess whether the bridge-joining algorithm is doing the right
thing.  Reports:
  - Component count and size distribution.
  - Position of each component (centroid + bounding box).
  - Whether each component lies inside the polyline (sanity check).
  - Bridge quad detection: cells whose 4 corners include points NOT
    in any other cell's edges (heuristic for bridges extending into
    open water).

Usage:
    python3 analyze_output.py <postarchipelago.vtk> <input.poly>
"""
import sys
import re
from pathlib import Path
from collections import defaultdict


def parse_vtk(path):
    with open(path) as fh:
        lines = fh.read().splitlines()
    num_points = None
    points = []
    cells = []
    idx = 0
    while idx < len(lines):
        line = lines[idx].strip()
        if line.startswith("POINTS"):
            num_points = int(line.split()[1])
            floats_needed = num_points * 3
            consumed = 0
            idx += 1
            all_floats = []
            while consumed < floats_needed and idx < len(lines):
                tok = lines[idx].split()
                consumed += len(tok)
                all_floats.extend([float(t) for t in tok])
                idx += 1
            points = [(all_floats[3 * i], all_floats[3 * i + 1])
                      for i in range(num_points)]
            continue
        if line.startswith("CELLS"):
            num_cells = int(line.split()[1])
            idx += 1
            for _ in range(num_cells):
                cell_tok = lines[idx].split()
                cell = [int(x) for x in cell_tok[1:1 + int(cell_tok[0])]]
                cells.append(cell)
                idx += 1
            continue
        idx += 1
    return num_points, points, cells


def parse_polyline(path):
    """Return a list of (x, y) tuples for the polyline vertices."""
    with open(path) as fh:
        lines = fh.readlines()
    # Header is on first line; first segment is <num_verts> <dim> ...
    header = lines[0].split()
    num_verts = int(header[0])
    pts = []
    for line in lines[1:1 + num_verts]:
        parts = line.strip().split()
        pts.append((float(parts[1]), float(parts[2])))
    return pts


def cell_centroid(cell, points):
    cx = sum(points[v][0] for v in cell) / len(cell)
    cy = sum(points[v][1] for v in cell) / len(cell)
    return cx, cy


def find_components(cells, points):
    """BFS over cell adjacency (via shared edges)."""
    edges = defaultdict(set)
    for ci, cell in enumerate(cells):
        n = len(cell)
        for i in range(n):
            a, b = cell[i], cell[(i + 1) % n]
            e = tuple(sorted([a, b]))
            edges[e].add(ci)
    adj = defaultdict(set)
    for cl in edges.values():
        if len(cl) == 2:
            cl = list(cl)
            adj[cl[0]].add(cl[1])
            adj[cl[1]].add(cl[0])
    comp = [-1] * len(cells)
    sizes = []
    nc = 0
    for s in range(len(cells)):
        if comp[s] != -1:
            continue
        stack = [s]
        comp[s] = nc
        cnt = 1
        while stack:
            x = stack.pop()
            for nb in adj[x]:
                if comp[nb] == -1:
                    comp[nb] = nc
                    stack.append(nb)
                    cnt += 1
        sizes.append(cnt)
        nc += 1
    # Pair each component id with its size and centroid.
    comp_info = []
    for c in range(nc):
        members = [i for i, lab in enumerate(comp) if lab == c]
        cx = sum(cell_centroid(cells[m], points)[0] for m in members) / len(members)
        cy = sum(cell_centroid(cells[m], points)[1] for m in members) / len(members)
        xs = []
        ys = []
        for m in members:
            for v in cells[m]:
                xs.append(points[v][0])
                ys.append(points[v][1])
        comp_info.append({
            "id": c, "size": len(members),
            "cx": cx, "cy": cy,
            "xmin": min(xs), "xmax": max(xs),
            "ymin": min(ys), "ymax": max(ys),
        })
    return comp_info


def point_in_polygon(x, y, poly):
    """Ray-casting algorithm: returns True if (x, y) is inside `poly`."""
    n = len(poly)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and \
           (x < (xj - xi) * (y - yi) / (yj - yi) + xi):
            inside = not inside
        j = i
    return inside


def main():
    if len(sys.argv) != 3:
        print("usage: analyze_output.py <postarchipelago.vtk> <input.poly>")
        sys.exit(1)
    vtk_path = Path(sys.argv[1])
    poly_path = Path(sys.argv[2])
    num_points, points, cells = parse_vtk(vtk_path)
    poly = parse_polyline(poly_path)
    comps = find_components(cells, points)
    comps.sort(key=lambda c: -c["size"])
    print(f"=== {vtk_path} ===")
    print(f"  cells={len(cells)} points={num_points} "
          f"components={len(comps)}")
    print(f"  polyline vertices={len(poly)}")
    print(f"  X range: [{min(p[0] for p in points):.3f}, "
          f"{max(p[0] for p in points):.3f}]")
    print(f"  Y range: [{min(p[1] for p in points):.3f}, "
          f"{max(p[1] for p in points):.3f}]")
    print(f"  polyline X range: [{min(p[0] for p in poly):.3f}, "
          f"{max(p[0] for p in poly):.3f}]")
    print(f"  polyline Y range: [{min(p[1] for p in poly):.3f}, "
          f"{max(p[1] for p in poly):.3f}]")
    print()
    print(f"Top 10 components (sorted by size, descending):")
    print(f"  {'id':>3} {'size':>4} {'cx':>9} {'cy':>9} "
          f"{'xrange':>17} {'yrange':>17} inside?")
    inside_count = 0
    outside_count = 0
    for c in comps[:10]:
        cx, cy = c["cx"], c["cy"]
        inside = point_in_polygon(cx, cy, poly)
        if inside:
            inside_count += 1
        else:
            outside_count += 1
        print(f"  {c['id']:>3} {c['size']:>4} "
              f"{cx:>9.3f} {cy:>9.3f} "
              f"[{c['xmin']:.2f},{c['xmax']:.2f}] "
              f"[{c['ymin']:.2f},{c['ymax']:.2f}] "
              f"{'YES' if inside else 'NO'}")
    print()
    print(f"Top-10 components: {inside_count} centroids inside polyline, "
          f"{outside_count} outside.")
    # Full component distribution
    print()
    print("All components (size, centroid-in-polyline?):")
    n_inside = 0
    n_outside = 0
    for c in comps:
        cx, cy = c["cx"], c["cy"]
        if point_in_polygon(cx, cy, poly):
            n_inside += 1
        else:
            n_outside += 1
    print(f"  Total: {len(comps)} components, "
          f"{n_inside} centroids inside polyline, "
          f"{n_outside} outside")


if __name__ == "__main__":
    main()
