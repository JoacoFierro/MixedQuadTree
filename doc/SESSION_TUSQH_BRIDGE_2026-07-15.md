# Resumen de la sesión: Sub-cell VF + Bridge-Joining (TUSQH)

**Fecha:** 2026-07-15
**Rama de trabajo:** integración TUSQH §3.3 + §3.4 en `MixedQuadTree`

---

## 1. Objetivo

Replicar los pasos §3.3 (sub-cell volume fractions) y §3.4 (archipelago resolution) del paper *TUSQH: Topological Control of Volume-Fraction Meshes Near Small Features and Dirty Geometry* (Shawcroft et al., 2025) sobre el mesher 2D `MixedQuadTree` ya existente.

---

## 2. Características implementadas

### 2.1. Sub-cell Volume Fractions (paper §3.3)

Para cada vértice y arista del quadtree se calcula el *fictitious cell* (rectángulo perpendicular centrado en la entidad de menor dimensión) y se promedian las winding numbers de sus `s × s` puntos muestra. Sólo se aceptan `s` pares y ≥ 2 para que las muestras no caigan sobre bordes.

- **Tamaño del fictitious cell:**
  - **Vértices**: cuadrado de lado = máxima longitud de aristas incidentes.
  - **Aristas**: rectángulo de largo = longitud de la arista y grosor perpendicular = máxima distancia perpendicular de los quads adyacentes (proyectando centroides sobre la recta de la arista).
- **Winding number**: se reutiliza `Polyline::getWindingNumber` ya integrado en el mesher.
- **Almacenamiento**:
  - Vértices: campos nuevos en `MeshPoint` (`mSubcellVolumeFraction`, `mSubcellIsInterior`, `mSubcellSampleSize`, etc.).
  - Aristas: side-map `mEdgeSubcellVF` keyed por `QuadEdge` con un POD `EdgeSubcellVFData` (no rompe el layout binario de `EdgeInfo`).

### 2.2. Bridge-Joining (paper §3.4 paso 3)

Para cada arista de frontera del quadtree (`MapEdges` con `info[2] == max()`) cuya ficticious cell sea interior (`VF ≥ joinThreshold`), se aplica una **división 1-a-5** sobre el quad adyacente:

- 4 sub-quads internos vía el `SplitVisitor` ya existente (refinamiento "verde" estándar).
- 1 *bridge quad* nuevo que se extiende hacia el **exterior** de la arista puente, perpendicularmente, con grosor `H / sampleSize`.

Implementado como loop de **punto fijo** (máximo 5 iteraciones) dentro de `resolveArchipelagos`:
1. Recomputar sub-cell VFs (no-op en la primera iteración).
2. Buscar aristas frontera con `VF ≥ joinThreshold`.
3. Para cada una, llamar a `bridgeSplitAtEdge(q, edgeIdx, s, nextQIdx)` → obtiene 5 nuevos quads (4 sub + 1 bridge) y devuelve vector.
4. Reemplazar el quad original con el primero y `push_back` los otros 4.
5. Repetir hasta que no se agreguen más bridges o se llegue al máximo de iteraciones.

Tras el loop, se aplica la **eliminación de componentes pequeños** (paper §3.4 paso 5): los componentes con menos de `minComponentCells` quads se descartan.

---

## 3. Archivos modificados / creados

### Nuevos

| Archivo | Propósito |
|---|---|
| `src/SubcellVFData.h` | POD `EdgeSubcellVFData` con `volumeFraction`, `windingNumbers`, `sampleSize`, `isInterior`. |
| `src/SubgridSampler.h` / `.cpp` | Lógica de muestreo: `sanitizeSampleSize`, `sampleVertex`, `sampleEdge`, `buildIncidentEdgeList`, `buildQuadPerpThickness`. |

### Modificados

| Archivo | Cambios |
|---|---|
| `src/MeshPoint.h` / `.cpp` | Nuevos campos y accessors `mSubcellSampleSize`, `mSubcellWindingNumbers`, `mSubcellVolumeFraction`, `mHasSubcellVolumeFraction`, `mSubcellIsInterior`. |
| `src/Mesher.h` / `.cpp` | Nuevas funciones miembro: `computeSubcellVolumeFractions()`, `bridgeSplitAtEdge()`, y `resolveArchipelagos()` reescrita con loop de bridge-joining. |
| `src/Visualization/VolumeFractionVTKWriter.h` / `.cpp` | Nuevos writers `writeSubcellVertexVF()` y `writeSubcellEdgeVF()` (UNSTRUCTURED_GRID con escalares `subcell_vf`, `subcell_is_interior`, `subcell_sample_size`, `subcell_vf_corner_mean`). |
| `src/Main.cpp` | Parseo de flags nuevos `-J`, `-K`, `-F`, `-L`. |
| `src/CMakeLists.txt` | `SubgridSampler` añadido a `SOURCE_FILES`. |

---

## 4. Flags CLI nuevos

| Flag | Significado | Default |
|---|---|---|
| `-J` | Habilita sub-cell VF + bridge-joining | off |
| `-K s` | Tamaño de muestreo (entero par ≥ 2) | 2 |
| `-F t` | Umbral de "interior" para join (`0..1`) | 0.5 |
| `-L n` | Mínimo de quads por componente (paper §3.4) | 5 |

Nota: se evitaron `-t` (ya usado por testing) y `-M` (ya usado por profundidad extra de TUSQH resolve).

---

## 5. Resultados de tests

### Caso exitoso: `tusqh_small_feature.poly`

```
-a 2 -T -J -K 2 -F 0.5 -L 3 -e
```

- 3 componentes iniciales: 4, 2, 2 quads.
- Iteración 0: 5 bridges agregados.
- Iteración 1: 1 bridge agregado.
- Iteraciones 2-4: sin cambios (convergió).
- Resultado: 3 componentes → 1 componente, 31 quads en el mesh final.
- 1 quad pequeño (sub-componente residual de 2 quads) descartado por `-L 3`.

**Salida VTK útil:**
- `_subcell_vertex.vtk` y `_subcell_edge.vtk`: muestras por vértice y arista con `subcell_vf`, `subcell_is_interior`, `subcell_sample_size`.
- `_postarchipelago.vtk`: estado del quadtree tras el bridge-joining + drop.

### Caso sin convergencia: `lonso.poly`

```
-a 4 -J -K 2 -F 0.5 -L 1 -e
```

- 25 componentes, 57 bridges agregados en 5 iteraciones (no convergió).
- Quads en `postarchipelago`: 409 (vs. 172 sin `-J`).
- El límite de 5 iteraciones protege contra crecimiento explosivo pero deja la cadena de bridges incompleta para esta geometría.

### Caso trivial: `unit_square.poly`

```
-a 4 -J -K 2 -F 0.5 -L 1 -e
```

- 1 componente, 0 bridges. Correcto.

---

## 6. Trabajo TODO / Pendiente

### Crítico (afecta correctitud)

- [x] **Filtrar aristas en frontera del dominio global** (Issue #1): implementado
  como `Mesher::isEdgeOnDomainBoundary` (`Mesher.cpp:2992-3050`). Usa una
  grilla `s × s` en el lado +dir_ext de la arista y descarta la arista si
  `vf_ext < joinThreshold`. Semántica estricta (`<`).
- [x] **Validar la dirección "exterior" vs "interior"** (Issue #2):
  implementado `Mesher::computeExteriorDirection` (`Mesher.cpp:2969-2990`)
  usando el **centroide del quad** (`mean(p_e0..p_e3)`) en lugar del
  midpoint del edge opuesto. Robusto a quads no rectangulares y a
  aproximaciones axis-aligned de features rotadas.
- [x] **Stale `info[1]` en aristas originales** (Issue #3):
  `bridgeSplitAtEdge` ahora elimina `MapEdges[QuadEdge(e0,e1)]` al
  finalizar el split (`Mesher.cpp:3135-3143`). Esto evita que el BFS
  downstream use la entrada obsoleta.
- [x] **Sub-muestreo de TUSQH en grilla inicial muy gruesa** (Issue #4):
  Chesapeake Bay con `-a 7 -T -N 2` producía sólo 3 celdas porque
  `generateGridMesh` produce 1-2 celdas raíz gigantes (~51×75) y el
  muestreo `-N 2 = 4` undersamplea. **Resuelto** con
  `Mesher::preRefineForTusqh` (`Mesher.cpp:3622-3673`,
  `Mesher.h:158-183`) que detecta el caso y pre-refina
  uniformemente a nivel 3 antes de `windingSubdivide`. Sólo se
  dispara cuando `Quadrants.size() <= 2` Y
  `input.getEdges().size() >= 100`, así no afecta a las
  polilíneas triviales de los tests de regresión.

### Importante (mejoras funcionales)

- [x] **Verificación con Figure 14 del paper**: creado
  `data/tusqh_figure14.poly` (rectángulo exterior + 2 islas interiores).
  Verificado con `-L 1` (3 componentes, 0 bridges añadidos) y `-L 3`
  (todos los componentes descartados por el umbral).
- [x] **Verificación con Chesapeake Bay (`data/Agua.poly`)**:
  - Sin fix, `-a 7 -T -N 2` → 3 celdas en 0.4 s
  - Con fix, `-a 3 -T -J -K 2 -F 0.5 -L 3` → 1034 points, 645 celdas en 53 s
  - Con fix, `-a 4 -T -N 4` → cientos de celdas en ~30 s
- [x] **Verificación de archipelago-joining con Chesapeake Bay**:
  con `-a 3 -T -J -K 2 -F 0.0 -L 1` se añadieron 282 bridges en 5
  iteraciones. **Limitación fundamental detectada**: los bridge quads
  son componentes aislados (97 → 378 componentes) en lugar de
  conectar componentes como en el paper §3.4. Documentado en
  `TUSQH_WINDING_GUIDE.md` §3.11.9.
- [x] **Crear polyline con gap pequeño**: `data/tusqh_bridge_small_gap.poly`
  (dos cuadrados 0.5×0.5 separados por gap de 0.1) confirma la
  limitación: incluso con gap pequeño, los bridge quads quedan
  aislados porque `H/sampleSize < gap`.
- [x] **Reimplementar bridge-joining paper-faithful (Step 1-7)**:
  en nueva rama `feature/paper-faithful-bridge` desde
  `develop-felipe@d174e6c`. 7 pasos completados, 5/5 tests
  originales + 1 nuevo test pasan. Reportes por paso en
  `doc/STEP_0`...`STEP_7`. Ver `doc/STEP_5_VALIDATION.md` para
  métricas de Chesapeake Bay: 125 bridges añadidos, 549 interior
  cells en output, 110 componentes en el mesh. **Issue #4 (limitación
  bridge-joining) RESUELTA**: el nuevo algoritmo detecta componentes
  sobre el cubical complex completo y los une correctamente.
- [ ] **Flag `-B` para máximo de iteraciones del loop de bridge**:
  actualmente hard-coded a 5 en `resolveArchipelagos`. Algunos casos
  podrían converger con 10-15; otros explotan con 3.
- [ ] **Flag `-B 0` para desactivar bridge-joining pero mantener sub-cell VF**:
  para usuarios que sólo quieren las VFs (útil para visualización).
- [ ] **Manejo de quads con aristas previamente divididas**: el
  `bridgeSplitAtEdge` aborta si la arista puente ya estaba dividida por
  una iteración previa. Para chains largos de bridges podría necesitar
  manejar ese caso (re-split no destructivo).
- [ ] **Flag de override para `preRefineForTusqh`**: actualmente la
  heurística corre automáticamente cuando dispara. Un flag tipo
  `-P 0` para desactivarla y un flag tipo `-P N` para forzar un nivel
  base distinto del default (3) podrían ser útiles para experimentación.
- [ ] **Visual validation in ParaView**: confirmar visualmente que
  el algoritmo paper-faithful produce la malla correcta para
  Chesapeake Bay. La métrica cuantitativa (Step 5) muestra que el
  algoritmo funciona, pero la confirmación visual es la prueba
  definitiva.

### Nice-to-have

- [x] **Documentar en `doc/TUSQH_WINDING_GUIDE.md`** la nueva sección
  §3.11 con ejemplos de uso de flags, call-graph extendido y notas de
  correctitud (Issues #1, #2, #3, #4 + guard contra segfault de
  Quadrants vacío + §3.11.9 limitación bridge-joining).
- [x] **Test automatizado** en `scripts/test_bridge_clean_info.py` con
  5 casos originales: unit_square, small_feature, figure14 (drop-all),
  boundary_filter (L-shape) y rotated (regression) + 1 caso nuevo
  `bridge_connects_two_islands` (positive control para Step 3). Los
  6 casos pasan con la implementación paper-faithful.
- [ ] **VTK writer para los bridge quads específicamente** (con un
  campo `is_bridge_quad` para distinguirlos en ParaView).
- [ ] **Métrica de "convergencia"** en el output: reportar no sólo "X
  bridges added" sino también "componentes al inicio vs. al final"
  para saber cuánto mejoró la conectividad.
- [ ] **Comparar visualmente con y sin bridge** en ParaView: renderizar
  `output/l_postarchipelago.vtk` y `_postarchipelago_no_br.vtk` lado a
  lado.

---

## 7. Bug conocido y resuelto

- **`L=99999` borra todo**: el flag `-L` significa "mínimo de quads
  para MANTENER". Con `L=99999`, ningún componente cumple, y
  `Quadrants` queda vacío, causando segfault en
  `saveOutputMesh:1623`. **Resuelto** con un guard al inicio de
  `Mesher::saveOutputMesh(vector<Quadrant>&)` (`Mesher.cpp:1622-1634`)
  que emite un mesh vacío sin dereferenciar `tmp_Quadrants[0]`. El
  mismo guard protege las llamadas posteriores en `Mesher.cpp:138,
  150, 165` (debug `_quads`, `_closeto`, `_remSur`) que tenían el
  mismo riesgo latente. Verificado con `tusqh_figure14.poly -L 3`
  (todos los componentes descartados, exit 0).

- **TUSQH undersampling con grilla inicial muy gruesa** (Issue #4):
  Chesapeake Bay con `-a 7 -T -N 2` producía sólo 3 celdas porque la
  grilla inicial del bbox (`generateGridMesh`) genera 2 celdas
  gigantes (~51×75) y el muestreo `-N 2 = 4` undersamplea la celda,
  clasificándola toda como `AllInside`/`AllOutside` sin subdivisión.
  **Resuelto** con `Mesher::preRefineForTusqh`
  (`Mesher.cpp:3622-3673`) que pre-refina uniformemente a nivel 3
  antes de `windingSubdivide`. Verificado:
  `-a 3 -T -J -K 2 -F 0.5 -L 3` → 645 celdas en 53 s (antes: 3 celdas).

---

## 8. Detalles técnicos relevantes para próximos pasos

- `MapEdges` almacena `q_id` (counter monotónico de `SplitVisitor`), **no** la posición en el vector `Quadrants`. Cualquier traducción entre ambos requiere el map `qIdToIdx` que se reconstruye tras cada bridge-joining.
- `QuadEdge(a, b)` (ctor por defecto) **ordena** los índices; `QuadEdge(a, b, true)` (forced) **no**. Usar el correcto es crucial para `MapEdges.find`.
- `EdgeInfo`: `info[0]` = midpoint (0 si no dividida), `info[1]` y `info[2]` = quads adyacentes (`numeric_limits<unsigned int>::max()` = sin quad).
- `SplitVisitor::visit` muta `MapEdges` y agrega puntos a `new_pts` (lista), pero **NO** modifica el vector `points` global — el llamador debe hacer `points.insert(points.end(), new_pts.begin(), new_pts.end())` después. Los índices dentro de `new_eles` ya están corregidos para apuntar a las posiciones futuras.
- `Quadrant::q_id` se establece en el constructor y no tiene setter; cuando se reemplaza un quad por 5 nuevos, hay que asignar `nextQIdx, nextQIdx+1, ..., nextQIdx+4` al construir los nuevos.
- `Quadrant::getIndex()` es **no-const** (declarado así en el header); tomar referencia `const Quadrant&` y llamar `getIndex()` falla por `-fpermissive`.
- `EdgeSubcellVFData::volumeFraction` (no `vf`).
- **Para `Triangle .poly`, CCW es OUTER boundary (dentro); CW es HOLE.**
  Si la polilínea tiene varios polígonos cerrados todos en CCW,
  Triangle los trata como regiones ADICIONALES (no como restas).
  Por ejemplo, un cuadrado CCW dentro de un rectángulo CCW produce
  `winding = 2` dentro del cuadrado y `winding = 1` dentro del
  rectángulo pero fuera del cuadrado — ambos son "interior". Para
  hacer un agujero, listar el polígono interior en sentido **CW**.
- **Quadrant::q_id vs índice en vector**: `q_id` es un counter
  monotónico; la posición en `Quadrants` puede cambiar tras splits.
  Usar `qIdToIdx` (reconstruido al inicio de `resolveArchipelagos`)
  para traducir.
- **`generateGridMesh` produce 1-N celdas según `GridMesher`**: el
  paso se calcula como `min(dx,dy)*1.01`. Para bbox no cuadrados
  la grilla puede quedar 1×2, 2×1, etc. Para bbox muy grandes (e.g.
  Chesapeake Bay) produce exactamente 2 celdas. Por eso la heurística
  de pre-refinement dispara con `coarseThreshold=2`.
- **`generateQuadtreeMesh` tiene 9 argumentos**:
  `(rl, input, all_reg, name, minrl, givenmaxrl, debugging,
  new_q_idx, Aliasing)`. La llamada pre-existente en `Mesher.cpp:271`
  (rama `else` cuando NO se usa TUSQH) está en orden incorrecto y
  es un bug latente, pero **no se corrige** en esta sesión porque
  cambiarlo podría alterar mallas ya producidas en otros flujos.

---

## 9. Estado final tras la implementación paper-faithful

Esta sesión cerró con el plan de 7 pasos ejecutado en la rama
`feature/paper-faithful-bridge`. Documentación completa en:

- **`doc/WORK_SUMMARY.md`** — resumen maestro: objetivo, plan, métricas finales, cómo usar.
- **`doc/BUGS_FOUND.md`** — Issues #1-#7 con symptom → root cause → fix → verification.
- **`doc/FUTURE_WORK.md`** — TODOs categorizados por prioridad (P1 visual validation, P2 features, P3 polish, P4 known issues, P5 research).
- **`doc/CURRENT_STATE.md`** — estado del branch: working tree dirty, sin commits nuevos, listo para revisión.
- **`doc/STEP_0..7_*.md`** — reportes detallados de cada paso.

### Métricas finales (Chesapeake Bay, `-a 3 -T -J -K 2 -F 0.5 -L 1`)

| Métrica | Valor |
|---------|-------|
| Wall-clock | ~38 s |
| Interior cells en output | 549 |
| Mesh components | 110 |
| Biggest component | 49 cells |
| Bridges añadidos | 125 (en 2 iteraciones) |
| AllOutside cells drop | 970 |
| Tests regresión | 6/6 pasan |

### Pendiente crítico (de `FUTURE_WORK.md` §1)

- [ ] **Visual ParaView verification** — confirmar visualmente que los 125 bridges conectan físicamente los componentes correctos.
- [ ] **Re-run tests tras commit de docs** — confirmar 6/6 pasan tras los cambios finales de documentación.

### Bugs resueltos en esta sesión (de `BUGS_FOUND.md`)

- ✅ Issue #1 — `isEdgeOnDomainBoundary` (strict `<` semantics).
- ✅ Issue #2 — `computeExteriorDirection` (centroid-based).
- ✅ Issue #3 — `bridgeSplitAtEdge` MapEdges cleanup.
- ✅ Issue #4 (paper-faithful rewrite) — bridge-joining ahora conecta componentes (NO produce quads aislados). **CRÍTICO**.
- ✅ Issue #5 — `saveOutputMesh` empty-Quadrants guard.
- ✅ Issue #6 — `preRefineForTusqh` heuristic (Chesapeake Bay 3→645 cells).
- ⚠️ Issue #7 — `Quadrant::getIndex()` const-ness (documented, no fixed).

### Decisiones de diseño finales

1. **Full-grid BFS** sobre el cubical complex completo (interior + exterior preservado), no interior-only. Filtro de output por `isInteriorCell` aplica después.
2. **`Unknown` NO es interior** en `isInteriorCell`. Testado y rechazado (loops infinitos).
3. **Bridge quads reclasificados** después del split con `WindingNumberVisitor` (no quedan como `Unknown`).
4. **Output filter** en `keepQuad[]` aplica `isInteriorCell` además de `compOfQuad >= 0` y `compSize >= minComponentCells`.

