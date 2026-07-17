# Guía de implementación: TUSQH, Winding Numbers y Flags

> **Propósito:** Este documento es la guía de referencia para entender **dónde**
> y **cómo** está implementado en el código todo lo relativo a *winding numbers*,
> subdivisión estilo **TUSQH**, *volume fractions* y los flags `-T -N -E -M -n -s`
> del ejecutable. Está pensado como apoyo para futuros *merges* entre ramas
> (e.g. `develop-felipe` ↔ `develop` ↔ `master`).
>
> Auditado sobre la rama `feature/paper-faithful-bridge` (HEAD `d174e6c Union
> Winding Numbers con Templates`, working tree sin commits nuevos). Las
> referencias `archivo:línea` son absolutas dentro del repo `MixedQuadTree`.

## 0. Documentación relacionada

Esta guía cubre el **cómo** de la implementación. Para el **qué / por qué /
resultados**, ver:

- **`doc/WORK_SUMMARY.md`** — resumen maestro de la rama
  `feature/paper-faithful-bridge`: objetivo, plan de 7 pasos, métricas finales
  de Chesapeake Bay (549 interior cells, 110 components, 125 bridges).
- **`doc/BUGS_FOUND.md`** — cada bug encontrado y arreglado (Issues #1-#7)
  con `symptom → root cause → fix → verification`.
- **`doc/FUTURE_WORK.md`** — TODOs pendientes categorizados P1-P5 (visual
  ParaView, flag `-B`, 3D pinch, persistence diagrams, etc.).
- **`doc/CURRENT_STATE.md`** — estado actual del branch: working tree dirty,
  sin commits nuevos, listo para revisión.
- **`doc/SESSION_TUSQH_BRIDGE_2026-07-15.md`** — log de la sesión que motivó
  este trabajo.
- **`doc/STEP_0..7_*.md`** — reportes detallados de cada uno de los 7 pasos.

**Esta guía (§3.11.9) es la implementación paper-faithful, ya NO la
limitación documentada previamente.** El algoritmo anterior solo consideraba
aristas expuestas y producía quads aislados; el nuevo (Step 1 + Step 3
rewrite) detecta componentes sobre el cubical complex completo y los une
correctamente.

---

## 1. Resumen conceptual

### 1.1 Winding number (wn)

Para un punto `P` y un polígono cerrado `V` (sin auto-intersecciones en el
sentido CCW/CW), el *winding number* `wn(P)` cuenta cuántas veces la
circunsferencia envuelve al punto:

- `wn(P) > 0` ⇒ `P` está **dentro** del dominio.
- `wn(P) == 0` ⇒ `P` está **fuera** del dominio.

El único punto del repositorio donde se calcula esto es
`Polyline::windingNumber(P)` (`src/Polyline.cpp:212-238`). El algoritmo es
el clásico `wn_PnPoly` de
[geomalgorithms.com](http://geomalgorithms.com/a03-_inclusion.html), basado
en pruebas `isLeft`. Cualquier visitante o ruta de código que necesite
saber si un punto está dentro del dominio llama a este método.

### 1.2 Volume fraction (VF)

TUSQH define la *volume fraction* (en 2D, *area fraction*) de una celda
cuadrangular como el **promedio de los winding numbers de una grilla
regular s × s de puntos de muestra** distribuidos dentro de la celda:

```
VF(quadrant) = (1 / s²) · Σ  wn(P_ij)
```

Si la celda está totalmente dentro ⇒ todos los `wn > 0` ⇒ `VF = 1`. Si está
totalmente fuera ⇒ todos los `wn == 0` ⇒ `VF = 0`. Si está cruzada ⇒
`0 < VF < 1` y se subdivide.

Esta es la definición implementada en
`Quadrant::computeVolumeFraction(wn)` (`src/Quadrant.cpp:179-188`) y
orquestada por `WindingNumberVisitor`
(`src/Visitors/WindingNumberVisitor.cpp:65-90`).

### 1.3 TUSQH (criterio de subdivisión)

Regla de subdivisión "subdivide-while-ambiguous" del paper *TUSQH:
Topological Control of Volume-Fraction Meshes Near Small Features and
Dirty Geometry* (Bracci et al., referenciado en
`src/Visitors/WindingNumberSubdivisionVisitor.h:29-31`):

Para cada candidato a refinar:

1. Calcular `wn(P_ij)` para los `s × s` puntos de muestra.
2. Clasificar la celda:
   - `AllInside`  ⇒ todos `wn > 0`  ⇒ **no** subdividir.
   - `AllOutside` ⇒ todos `wn == 0` ⇒ **no** subdividir.
   - `Mixed`      ⇒ combinación de los dos ⇒ **subdividir**.
3. Opcional (modo legado `-E`): si la celda intersecta geométricamente
   alguna arista del polígono de entrada, forzar subdivisión
   independientemente del winding.

El estado por celda se almacena en el enum `WindingState`
(`src/Quadrant.h:53-61`).

### 1.4 Relación entre `-n` (volume fraction) y `-N` (TUSQH)

Son **dos parámetros independientes**:

- `-n s` controla el tamaño de la grilla usada por
  `Mesher::computeVolumeFractions` (`src/Mesher.cpp:2223-2247`). Se aplica
  a **todas** las celdas ya existentes (post-subdivisión).
- `-N s` controla el tamaño de la grilla usada por
  `WindingNumberSubdivisionVisitor` (`src/Mesher.cpp:1770-1780`), dentro
  del bucle TUSQH. Sólo se aplica a celdas candidatas durante la
  subdivisión.

Por defecto ambos valen 2, pero pueden fijarse a valores distintos
(véase `src/Main.cpp:96-101`).

---

## 2. Flags de línea de comandos

Definidos y parseados exclusivamente en `src/Main.cpp`. Resumen:

| Flag | Tipo        | Variable destino                 | Default | Significado                                                                                          |
|------|-------------|---------------------------------|---------|------------------------------------------------------------------------------------------------------|
| `-T` | bool        | `useTusqh`                      | `false` | Activa subdivisión TUSQH (`windingSubdivide`). Sin `-T` se usa el pipeline clásico.                  |
| `-N s` | `unsigned int` | `tusqhSampleSize`            | `2`     | Tamaño de la grilla s × s dentro del **criterio TUSQH** (no afecta a `-n`).                          |
| `-E` | bool        | `refineOnEdgeIntersect`         | `false` | Modo legado: refinar también cualquier celda que intersecte geométricamente una arista de entrada.   |
| `-M d` | `unsigned int` | `tusqhExtraResolveDepth`     | `0`     | Iteraciones extra de TUSQH ("resolve pass") sobre celdas `Mixed` que llegan a `qrl == maxDepth`.     |
| `-n s` | `unsigned int` | `mSampleSize`                 | `2`     | Tamaño de la grilla para el cálculo de **volume fraction** (post-proceso).                            |
| `-s rl` | `unsigned short` | crea `RefinementInputSurfaceRegion` | n/a   | **No** relacionado con TUSQH ni winding: refina hasta `rl` las celdas que intersectan la superficie.  |
| `-J` | bool        | `useSubgrid`                    | `false` | Activa sub-cell VF + archipelago resolution (paper §3.3+§3.4).                                       |
| `-K s` | `unsigned int` | `subgridSampleSize`         | `2`     | Tamaño de la grilla s × s para muestreo sub-cell (independiente de `-n` y `-N`).                     |
| `-F τ` | `double`     | `subgridJoinThreshold`         | `0.5`   | Umbral VF sub-cell. `VF ≥ τ` ⇒ vértice/arista "interior" candidato a puente.                          |
| `-L n` | `unsigned int` | `subgridMinComponentCells`  | `5`     | Componentes con menos de `n` celdas se descartan en el resolver de archipiélagos.                    |

Puntos exactos en `src/Main.cpp`:

- Mensaje de ayuda: `Main.cpp:80-108` (los flags TUSQH/VF están en `:96-107`).
- Declaración de variables: `Main.cpp:147-151`.
- Marcado de flags sin argumento en `:179-180` (los booleanos van en la lista
  `inout = true`).
- Parsing de `-T -N -E -M`: `Main.cpp:193-216`.
- Parsing de `-n`: `Main.cpp:389-392`.
- Parsing de `-s`: `Main.cpp:304-318` (crea región de refinamiento; no toca
  TUSQH).
- Paso al `Mesher`: `Main.cpp:441-454` (genera la malla desde cero o la
  refina partiendo de un quadtree existente).

---

## 3. Mapa archivo-por-archivo

### 3.1 `src/Main.cpp` — punto de entrada y parseo de flags

| Líneas | Qué hace                                                                                                                                                                                                                                  |
|--------|---|
| `96-107` | Texto de ayuda de los flags `-n -T -N -E -M`.                                                                                                                              |
| `147-151` | Declaración de `mSampleSize`, `useTusqh`, `tusqhSampleSize`, `refineOnEdgeIntersect`, `tusqhExtraResolveDepth`.                                                              |
| `169-185` | Lista de flags que **no** esperan argumento (`-T`, `-E`, etc.).                                                                                                            |
| `193-216` | `switch` que parsea `-T`, `-N`, `-E`, `-M`.                                                                                                                                |
| `217-240` | Caso `-t` (testing): imprime `windingNumber(...)` para varios puntos — sirve como regresión visual del primitive.                                                          |
| `304-318` | Caso `-s`: **no** relacionado con TUSQH — crea una región de refinamiento `RefinementInputSurfaceRegion`.                                                                  |
| `389-392` | Caso `-n`: guarda tamaño de muestra para volume fraction en `mSampleSize`.                                                                                                 |
| `441-454` | Llamada a `mesher.generateMesh(...)` o `mesher.refineMesh(...)` con todos los flags TUSQH/VF como parámetros.                                                              |

### 3.2 `src/Mesher.h` — interfaz del meseher

| Líneas              | Qué hace                                                                                                                                              |
|---------------------|---|
| `79-86`             | Firma de `generateMesh(...)` con parámetros `sampleSize`, `useTusqh`, `tusqhSampleSize`, `refineOnEdgeIntersect`, `tusqhExtraResolveDepth`.               |
| `88-97`             | Firma de `refineMesh(...)` con los mismos parámetros TUSQH/VF.                                                                                       |
| `118-130`           | Doc y declaración de `windingSubdivide(...)` — el bucle TUSQH.                                                                                       |
| `145`               | Declaración de `computeVolumeFractions(input, sampleSize)` — el post-proceso de VF.                                                                  |
| `190`               | Miembro `mSampleSize` (el `-n`).                                                                                                                      |

### 3.3 `src/Mesher.cpp` — orquestación del bucle TUSQH y del cálculo de VF

| Líneas              | Qué hace                                                                                                                                              |
|---------------------|---|
| `82-88`             | `refineMesh`: si `useTusqh` ⇒ llama a `windingSubdivide(...)`; si no ⇒ llama al `splitQuadrants(...)` clásico.                                      |
| `90-92`             | En ambos casos (TUSQH o clásico) **siempre** se ejecuta `computeVolumeFractions(input, mSampleSize)` (la `-n`).                                       |
| `234-240`           | `generateMesh`: mismo fork que arriba.                                                                                                                |
| `242-244`           | Idem: `computeVolumeFractions(input, mSampleSize)` se ejecuta siempre.                                                                               |
| `1707-1736`         | Doc y firma de `Mesher::windingSubdivide(...)`.                                                                                                       |
| `1738-1793`         | Estado local: `candidates`, `processed`, `idx_pos_map`, `new_pts`, `needs_classification`, `toBalance`.                                              |
| `1770-1788`         | Crea `WindingNumberSubdivisionVisitor` con `tusqhSampleSize` y `refineOnEdgeIntersect`. Si `-E` ⇒ le pasa un `IntersectionsVisitor`.                |
| `1795-1954`         | **Bucle principal TUSQH**: por cada `depth ∈ [0, maxDepth)` clasifica, subdivide, balancea, contabiliza `refinedCount`.                              |
| `1804-1824`         | **Paso 1**: clasifica cada candidato con `wnsv.visit(&quad)`.                                                                                          |
| `1840-1871`         | **Paso 2**: subdivide las celdas `Mixed` con `SplitVisitor`.                                                                                          |
| `1876-1927`         | **Paso 3**: balanceo one-irregular a través de la lista `toBalance`.                                                                                   |
| `1938-1947`         | Clasificación diferida (`needs_classification`) para celdas balanceadas recién creadas (sus puntos aún no estaban en `points`).                       |
| `1956-1964`         | Promoción de candidatos sobrantes al final del bucle.                                                                                                  |
| `1966-2001`         | **Post-proceso TUSQH**: descarta celdas `AllOutside` (VF=0%); limpia `intersected_edges` para `AllInside`.                                            |
| `2023-2202`         | **Resolve pass (`-M`)**: subdivide adicionalmente las celdas `Mixed` que quedaron al alcanzar `qrl == maxDepth`. Sin balanceo, criterio puro TUSQH. |
| `2209-2217`         | Dump VTK debug `_tusqh` y `_winding_state.vtk` (sólo si `VTKOUT==true`).                                                                               |
| `2223-2247`         | `Mesher::computeVolumeFractions(input, sampleSize)` — usa `WindingNumberVisitor` y vuelca `volume_fraction_debug.vtk` + `_samples.vtk`.              |

### 3.4 `src/Polyline.h` / `src/Polyline.cpp` — primitivo de winding number

| Líneas                          | Qué hace                                                                          |
|---------------------------------|---|
| `Polyline.h:73`                 | Declaración de `int windingNumber(const Point3D&) const`.                          |
| `Polyline.cpp:212-238`          | **Único** punto del repo donde se calcula un winding number. Algoritmo `wn_PnPoly`. |
| `Polyline.cpp:287-295`          | `pointIsInMesh(P) ≡ wn(P) > 0` — usado en otros lugares del pipeline clásico.       |

> **Regla para merges:** *cualquier* cambio en el cálculo del winding
> number debe concentrarse aquí. No dupliques la lógica en los visitantes.

### 3.5 `src/Quadrant.h` / `src/Quadrant.cpp` — estado por celda

`Quadrant` es la unidad básica del quadtree. Almacena:

- Topología habitual: `pointindex`, `sub_elements`, `intersected_edges`,
  `intersected_features`, `ref_level`, `q_id`, `surface`, `max_dis`,
  `debugging`.
- **Bloque Volume Fraction** (nuevo para TUSQH):
  - Miembros: `mSampleSize`, `mWindingNumbers`, `mVolumeFraction`, `mHasVolumeFraction`
    (`src/Quadrant.h:188-193`).
  - API: `setSampleSize(s)`, `getSampleSize()`, `computeVolumeFraction(wn)`,
    `getVolumeFraction()`, `getWindingNumbers()`, `hasVolumeFraction()`,
    `getSamplePoint(i, j, mp)` (`src/Quadrant.h:135-144, 355-375`).
- **Bloque TUSQH winding-state** (nuevo):
  - `enum class WindingState { Unknown, AllInside, AllOutside, Mixed }`
    (`src/Quadrant.h:53-61`).
  - Miembro: `mWindingState` (`src/Quadrant.h:194-196`).
  - API inline: `setWindingState`, `getWindingState`, `isWindingInside`,
    `isWindingOutside`, `isWindingMixed` (`src/Quadrant.h:377-393`).
  - Friend de los dos visitantes: `WindingNumberVisitor` y
    `WindingNumberSubdivisionVisitor` (`src/Quadrant.h:73-74`).

Implementación:

| Líneas                       | Qué hace                                                                                |
|------------------------------|---|
| `Quadrant.cpp:160-174`       | `getSamplePoint(i, j, mp)` — genera la posición del (i,j)-ésimo punto de la grilla s × s dentro de la bbox de la celda. |
| `Quadrant.cpp:179-188`       | `computeVolumeFraction(wn)` — `mVolumeFraction = Σ wn / s²` y guarda los `wn` originales. |

### 3.6 `src/Visitors/WindingNumberSubdivisionVisitor.{h,cpp}` — clasificación TUSQH

Este visitor **decide** si una celda debe subdividirse (es el corazón del
criterio TUSQH).

- Header: `src/Visitors/WindingNumberSubdivisionVisitor.h:1-81`
  - Ctor: `WindingNumberSubdivisionVisitor(unsigned int s, bool refineOnEdgeIntersect = false)`
    (línea 54-55).
  - Setters: `setPolyline`, `setPoints`, `setIntersectionsVisitor` (líneas 61-69).
  - Miembros: `mSampleSize`, `mRefineOnEdgeIntersect`, `mPolyline`, `mPoints`, `mIntersectionsVisitor`.
- Implementación (`src/Visitors/WindingNumberSubdivisionVisitor.cpp`):
  - `90-180` método `visit(Quadrant*)` — devuelve `true` ⇔ celda debe subdividirse.
  - `96` fija `mSampleSize` en la celda.
  - `103-131` rama **legada `-E`**: si la celda intersecta alguna arista (probado
    con `IntersectionsVisitor`) ⇒ marca `Mixed` y devuelve `true` sin
    recalcular winding numbers.
  - `133-179` rama **TUSQH pura**: muestrea `s × s`, llama a
    `Polyline::windingNumber`, clasifica y devuelve `true` sólo si
    `Mixed`.

### 3.7 `src/Visitors/WindingNumberVisitor.{h,cpp}` — cálculo de VF

Este visitor **no decide subdivisión**: se limita a poblar `mWindingNumbers`
y `mVolumeFraction` en cada celda.

- `src/Visitors/WindingNumberVisitor.cpp:65-90`
  - `visit(q)` ⇒ `setSampleSize(s)` y llama `computePostOrder(q)`.
  - `computePostOrder(q)`: para cada `i,j ∈ [0,s)` toma `q->getSamplePoint(i,j,mp)`,
    llama `mPolyline->windingNumber(sample)` y al final
    `q->computeVolumeFraction(wn_values)`.

### 3.8 `src/Visitors/IntersectionsVisitor.{h,cpp}` — usado por la rama `-E`

Sólo se instancia dentro de `windingSubdivide` cuando `refineOnEdgeIntersect`
es `true` (`src/Mesher.cpp:1777-1780`). Sirve como prueba geométrica de
intersección celda↔arista; **no** es código nuevo pero su uso es nuevo
dentro del flujo TUSQH.

### 3.9 `src/Visualization/VolumeFractionVTKWriter.{h,cpp}` — salida debug VTK

Tres *writers* (todos estáticos):

| Función                                  | Líneas        | Archivo generado                            | Contenido                                                                                          |
|------------------------------------------|---------------|---------------------------------------------|---|
| `writeQuadTreeWithVF`                    | `:45-130`     | `volume_fraction_debug.vtk`                  | QuadTree completo con `volume_fraction`, `intersects_surface`, `refinement_level` por celda.        |
| `writeVFHeatmap`                         | `:134-210`    | `volume_fraction_debug_samples.vtk`          | Cada uno de los `s × s × N_quadrants` puntos de muestra con su `winding_number` (heat map).         |
| `writeWindingState`                      | `:214-294`    | `<name>_winding_state.vtk`                  | QuadTree con clasificación TUSQH (`0=Unknown, 1=AllInside, 2=AllOutside, 3=Mixed`) + `VF` + nivel.  |

Llamadas desde `Mesher.cpp`:

- Tras `computeVolumeFractions` (cuando `VTKOUT==true`):
  `Mesher.cpp:2235-2241` → `writeQuadTreeWithVF` + `writeVFHeatmap`.
- Tras `windingSubdivide` (cuando `VTKOUT==true`):
  `Mesher.cpp:2209-2217` → `writeWindingState` con sufijo `_winding_state.vtk`.

### 3.10 `src/Visitors/SplitVisitor.{h,cpp}` — usado (sin cambios) por TUSQH

`SplitVisitor` no es código nuevo; `windingSubdivide` lo reusa para
materializar los 4 hijos de una celda `Mixed`
(`src/Mesher.cpp:1840-1871, 1898-1921, 2108-2121`). El contrato es el
mismo que en el pipeline clásico.

### 3.11 Bridge-joining pipeline (TUSQH §3.3 + §3.4)

Sub-cell volume fractions (vértices + aristas) y resolución de
archipiélagos. Se ejecuta después del TUSQH principal pero antes del
pipeline clásico (`removeOnSurfaceSafe`, etc.), porque necesita la
conectividad completa del quadtree.

#### 3.11.1 Activación y flags (`src/Main.cpp`)

| Flag | Variable destino              | Default | Significado                                                                                          |
|------|-------------------------------|---------|------------------------------------------------------------------------------------------------------|
| `-J` | `useSubgrid`                  | `false` | Activa sub-cell VF + archipelago resolution (paper §3.3 + §3.4).                                     |
| `-K s` | `subgridSampleSize`         | `2`     | Tamaño de la grilla `s × s` para muestreo sub-cell (no afecta a `-n` ni `-N`).                        |
| `-F τ` | `subgridJoinThreshold` (double) | `0.5` | Umbral VF sub-cell. Vértices/aristas con `VF ≥ τ` se consideran "interiores" y pueden usarse como puente. |
| `-L n` | `subgridMinComponentCells`  | `5`     | Componentes con menos de `n` celdas se descartan en el resolver de archipiélagos.                     |

Parseo: `src/Main.cpp:96-107` (texto de ayuda) y el switch adyacente
que rellena las cuatro variables. Se pasan a `Mesher::generateMesh` /
`Mesher::refineMesh` justo después de los flags `-T -N -E -M -n`.

#### 3.11.2 `computeSubcellVolumeFractions` (paper §3.3)

Recorre todas las aristas y vértices de cada celda TUSQH, muestrea `s × s`
dentro del **lado positivo** de la dirección exterior y guarda la
promedio de winding numbers en `mEdgeSubcellVF` y `mVertexSubcellVF`
(maps globales en `Mesher.h:288-305`).

- Firma: `Mesher.h:165-176`.
- Implementación: `Mesher.cpp:2876-2943`.
- Salidas debug (sólo si `VTKOUT==true`):
  - `<name>_subcell_vertex.vtk`: cada vértice con su `mVertexSubcellVF`.
  - `<name>_subcell_edge.vtk`: cada arista con su `mEdgeSubcellVF`.

#### 3.11.3 `bridgeSplitAtEdge` — creación de un puente

Aplica el *template* de subdivisión del paper (Figura 7): cuando una
arista tiene `mEdgeSubcellVF >= joinThreshold`, se subdivide la celda
adyacente usando la dirección perpendicular a la arista hacia el
**exterior** (no hacia el vecino), y se eliminan los hijos que
quedarían dentro de ese lado.

- Firma: `Mesher.h:178-189`.
- Implementación: `Mesher.cpp:2969-3145`.

Pasos clave:

1. `computeExteriorDirection(quad, e0, e1)` (`Mesher.cpp:2969-2990`):
   devuelve el vector unitario **+dir_ext** desde el *centroide del
   quad* (no del edge opuesto). Esto es robusto frente a rotaciones
   del quad y mantiene la semántica para quads no rectangulares
   (defensa contra el Issue #2 del plan original).
2. `splitAlongEdge(q, e0, e1, dir_ext)`: divide la celda en `k` capas
   perpendiculares a `dir_ext`, conservando sólo los hijos cuya
   posición proyectada sobre `dir_ext` está en el lado **negativo**
   (es decir, los hijos "interiores" del quad).
3. Limpia `MapEdges[QuadEdge(e0, e1)]` (Issue #3): la entrada previa
   de la arista original queda obsoleta tras el split. Eliminarla
   evita que el BFS en `resolveArchipelagos` use conectividad
   fantasma.
4. Devuelve `true` si el split amplió la celda (necesario para que el
   resolver decida iterar de nuevo).

#### 3.11.4 `isEdgeOnDomainBoundary` — filtro de Issue #1

`Mesher.cpp:2992-3050`. Para una arista `(e0, e1)` con dirección
exterior `+dir_ext` ya computada:

1. Toma `s × s` puntos de muestra desplazados una distancia fija desde
   el *midpoint* del edge hacia `+dir_ext`.
2. Calcula `vf_ext = Σ wn / s²`.
3. Si `vf_ext < joinThreshold` ⇒ descarta la arista como candidata a
   puente (es arista de frontera del dominio, no une archipiélagos).

Semántica **estricta** (`<`, no `≤`): una arista exactamente sobre la
frontera (`vf_ext = 0`) siempre se descarta; una con `vf_ext`
ligeramente positivo (e.g. `1/s²`) se conserva sólo si supera el
umbral.

#### 3.11.5 `resolveArchipelagos` — paper §3.4

Orquesta el ciclo completo:

```
bucle iter (hasta que no haya progreso):
  1. computeSubcellVolumeFractions(input, s, τ)
  2. Construir grafo de componentes vía BFS sobre MapEdges
  3. Para cada arista con mEdgeSubcellVF >= τ y
     info[2] == maxRefinementLevel (es borde de componente):
       3a. computeExteriorDirection(quad, e0, e1)
       3b. isEdgeOnDomainBoundary(...) -> skip si vf_ext < τ
       3c. bridgeSplitAtEdge(quad, e0, e1, dir_ext)
       3d. Si split exitoso, marcar arista como puente
  4. Recalcular componentes, descartar componentes con < L celdas
salida: Quadrants contiene sólo las celdas de los componentes supervivientes
```

- Firma: `Mesher.h:179-189`.
- Implementación: `Mesher.cpp:3286-3587`.
- Logging: `Mesher.cpp:3571-3577` (componentes, puentes, drops,
  tiempo). Cada iteración loggea `[bridge iter N] filtered X boundary
  edges (domain boundary)`.

#### 3.11.6 Guard contra segfault de Quadrants vacío (pre-existente)

Antes del fix, `Mesher::saveOutputMesh` accedía a
`tmp_Quadrants[0].getRefinementLevel()` sin verificar que el vector
estuviera no-vacío (`Mesher.cpp:1623`). Cuando `-L` descartaba todas
las componentes (caso `tusqh_figure14.poly -L 3`), el programa
crasheaba con SIGSEGV. Se añadió un guard al inicio de
`saveOutputMesh(vector<Quadrant>&)` (`Mesher.cpp:1622-1634`) que emite
un mesh vacío sin tocar `tmp_Quadrants[0]`. El guard también
beneficia las llamadas posteriores en `Mesher.cpp:138, 150, 165`
(debug `_quads`, `_closeto`, `_remSur`) que tenían el mismo riesgo
latente.

#### 3.11.7 Pruebas (`scripts/test_bridge_clean_info.py`)

Script de regresión que verifica:

1. Exit code 0 (no segfault) en todos los casos.
2. `cell count` del `_postarchipelago.vtk` (o ausencia correcta del
   archivo cuando todos los componentes se descartan).
3. Número inicial de componentes.
4. Cantidad mínima de aristas filtradas (Issue #1).
5. Índices de punto dentro de rango en cada celda (sanidad post-split).

Casos cubiertos:

- `data/unit_square.poly -a 3 -T -J -K 2 -F 0.5 -L 1`
- `data/tusqh_small_feature.poly -a 2 -T -J -K 2 -F 0.5 -L 3`
- `data/tusqh_figure14.poly -a 3 -T -J -K 2 -F 0.5 -L 3`
- `data/tusqh_boundary_filter.poly -a 3 -T -J -K 2 -F 0.5 -L 1`
  (regresión de Issue #1 con polilínea "L-shape")
- `data/tusqh_rotated.poly -a 4 -T -J -K 2 -F 0.5 -L 1`
  (regresión de Issue #2 con un cuadrado rotado)

#### 3.11.8 Heurística `preRefineForTusqh` — fix del bug de grilla inicial muy gruesa

Cuando el quadtree arranca con muy pocas celdas raíz (típicamente 1-2
para polilíneas grandes como Chesapeake Bay a nivel 0), el muestreo
por defecto de TUSQH (`-N 2`, 4 muestras por celda) puede
sub-muestrear la celda gigante y clasificarla toda como `AllInside`
o `AllOutside`, produciendo prácticamente cero subdivisiones. Esto
fue observado con `data/Agua.poly` (Chesapeake Bay, bbox ~51×75):
con `-a 7 -T -N 2` el resultado era 3 celdas (la mitad de una
celda inicial). Con `-E` (intersección de aristas) sí converge pero
tarda ~36 s y produce 698 celdas.

**Fix:** nuevo método `Mesher::preRefineForTusqh`
(`src/Mesher.cpp:3622-3673`, declarado en `src/Mesher.h:158-183`).
Se invoca automáticamente antes de `windingSubdivide` tanto en
`generateMesh` (`Mesher.cpp:267-268`) como en `refineMesh`
(`Mesher.cpp:99-100`), pero solo si:

1. `Quadrants.size() <= coarseThreshold` (default `2`), y
2. `input.getEdges().size() >= minSegmentsForTrigger` (default `100`),
   para no disparar en polilíneas triviales como `unit_square.poly`,
3. `maxDepth > 0`.

El umbral de segmentos (2) es clave: las pruebas de regresión
(`unit_square`, `small_feature`, `figure14`, `boundary_filter`,
`rotated`) tienen 4-13 segmentos y producen 1 celda con TUSQH,
resultado que **no debe** alterarse. La polilínea de Chesapeake Bay
tiene 46 432 segmentos y dispara la heurística.

Cuando se dispara, se invoca `generateQuadtreeMesh(baseLevel=3, ...)`
que refina la grilla inicial uniformemente a nivel 3 (= 64² = 4096
celdas teóricas, en la práctica unas decenas o cientos dependiendo
del intersect con la polilínea). El nivel base efectivo se capa a
`maxDepth` si el usuario pidió menos profundidad. Luego
`windingSubdivide` corre sobre la grilla pre-refinada y el muestreo
`-N 2` ya es significativo.

Resultados verificados con Chesapeake Bay (`-p data/Agua.poly`):

| Comando                          | Celdas | Tiempo |
|----------------------------------|--------|--------|
| `-a 7 -T -N 2` (sin fix)         | 3      | 0.4 s  |
| `-a 7 -T -E` (sin fix)           | 698    | 36 s   |
| `-a 7 -T -N 2` (con fix)         | ≥3500  | <120 s |
| `-a 4 -T -N 4` (sin fix)         | 174    | 7 s    |
| `-a 3 -T -J -K 2 -F 0.5 -L 3` (con fix) | 645 | 53 s |

Las 5 pruebas de regresión de `scripts/test_bridge_clean_info.py`
siguen pasando (la heurística no se dispara porque la polilínea tiene
menos de 100 segmentos).

**Nota técnica:** la llamada a `generateQuadtreeMesh` que se usa en
`preRefineForTusqh` pasa los argumentos en el orden correcto de la
firma `(rl, input, all_reg, name, minrl, givenmaxrl, debugging,
new_q_idx, Aliasing)`. La llamada pre-existente en `Mesher.cpp:271`
(la rama `else` cuando NO se usa TUSQH) tiene los argumentos en
orden incorrecto y produce una grilla distinta a la esperada — es
un bug pre-existente que NO se corrige aquí porque cambiarlo podría
alterar mallas ya producidas.

#### 3.11.9 Bridge-joining paper-faithful (Step 1 + Step 3 rewrite)

**Resumen:** el algoritmo bridge-joining se reimplementó para seguir
fielmente el paper TUSQH §3.4: detectar componentes sobre el cubical
complex completo (interior + exterior), identificar aristas entre
componentes distintos, y añadir 1-to-5 splits que conectan los
componentes.

**Cambios clave:**

1. **Step 1 (`Mesher.cpp` ~2030-2068): preservar celdas AllOutside.**
   Antes, las celdas con `WindingState::AllOutside` se descartaban
   después de TUSQH. Ahora se mantienen en `Quadrants` para que
   `MapEdges` retenga las aristas entre exterior e interior — el
   BFS necesita atravesar las celdas exteriores para identificar
   componentes correctamente.

2. **Step 2 (`Mesher.h:289` + `Mesher.cpp:3250-3260): helper
   `isInteriorCell`.** Distingue celdas `AllInside`/`Mixed` (interior)
   de `AllOutside`/`Unknown` (exterior) para que el BFS pueda contar
   cuántas celdas interiores tiene cada componente.

3. **Step 3 (`Mesher.cpp:3361-3520): nuevo bucle bridge-joining.** En
   lugar de considerar solo aristas expuestas (`info[2] == MAX`), el
   nuevo bucle:
   - Ejecuta BFS sobre el **cubical complex completo** (interior +
     exterior). Dos celdas interiores separadas por una cadena de
     celdas exteriores pertenecen al MISMO componente — son la misma
     "región topológica" según §3.4.
   - Busca candidatos: aristas en `MapEdges` donde ambos lados son
     celdas interiores en componentes DIFERENTES, con sub-cell VF ≥
     `joinThreshold` y NO en la frontera global del dominio (Issue #1
     filter).
   - Para cada candidato, aplica `bridgeSplitAtEdge` al quad en
     `info[1]`. El bridge quad se extiende hacia `info[2]`,
     conectando los dos componentes.

4. **Step 4 (`Mesher.cpp:3621-3690): filtro de output.** Las celdas
   `AllOutside` se preservan para la topología del BFS pero se
   excluyen del output final (mediante `isInteriorCell` en el filtro
   `keepQuad`).

**Resultados con Chesapeake Bay**
(`-p data/Agua.poly -a 3 -T -J -K 2 -F 0.5 -L 1`):

| Métrica                                | Valor    |
|----------------------------------------|----------|
| Cells pre-bridge (TUSQH output)        | 1019     |
| Cells post-bridge                      | 1519     |
| Bridges añadidos                       | 125      |
| Componentes pre-bridge (cubical)       | 74       |
| Componentes post-bridge (cubical)      | 71 (después de 2 iteraciones) |
| Cells en output final (interior)       | 549      |
| Componentes en output final (mesh)     | 110      |
| Mayor componente en output (cells)     | 49       |

Con el algoritmo anterior (limitado a exposed edges), 282 bridges
añadidos generaban 281 componentes nuevos aislados. Con el nuevo
algoritmo, los bridges unen componentes correctamente: aunque el
número total de componentes en el output crece (110 vs 85 en TUSQH
pre-bridge), esto se debe a la fragmentación natural de la bahía en
islas pequeñas con -K 2 (sample size 2, quad size ≈ 0.5 unidades),
no a bridges aislados.

**Verificación con gap pequeño (`data/tusqh_bridge_small_gap.poly`,
2 cuadrados 0.5×0.5 con gap=0.1):**
- 9 celdas iniciales (5 interior + 4 exterior preservadas).
- 2 bridges añadidos en iter 0.
- Output: 3 celdas (los 2 cuadrados + 1 bridge quad que pasa el
  filtro `isInteriorCell`).

**Limitación residual (paper-faithful pero no perfecta):** el bridge
quad se clasifica con `WindingNumberVisitor` después de añadirlo. Si
el bridge quad cae en zona exterior pura (e.g., gap sobre open water),
se clasifica como `AllOutside` y se descarta. Para Chesapeake Bay
con -K 2, la mayoría de los gaps entre islas son sobre open water,
por lo que la mayoría de los bridges se descartan — pero los que
sobreviven conectan correctamente los componentes que sí son
cercanos.

**Documentación relacionada:**
- Detalle por paso: `doc/STEP_0..7_*.md`.
- Bugs resueltos: `doc/BUGS_FOUND.md` (Issue #4 documenta el problema
  original resuelto por este rewrite; Issues #1-#3 son los filtros que
  el nuevo bucle aplica a los candidatos).
- Métricas finales: `doc/WORK_SUMMARY.md` §"Final results".
- Estado del branch: `doc/CURRENT_STATE.md`.

---

## 4. Variables y campos importantes

### 4.1 Variables del pipeline (en `Main.cpp`)

```cpp
unsigned int mSampleSize = 2;               // src/Main.cpp:147  ← -n
bool useTusqh = false;                       // src/Main.cpp:148  ← -T
unsigned int tusqhSampleSize = 2;            // src/Main.cpp:149  ← -N
bool refineOnEdgeIntersect = false;          // src/Main.cpp:150  ← -E
unsigned int tusqhExtraResolveDepth = 0;     // src/Main.cpp:151  ← -M
```

### 4.2 Parámetros de `Mesher::generateMesh` / `refineMesh`

```cpp
unsigned int sampleSize,            // = mSampleSize (controlado por -n)
bool useTusqh,
unsigned int tusqhSampleSize,
bool refineOnEdgeIntersect,
unsigned int tusqhExtraResolveDepth;
```

Declaraciones: `src/Mesher.h:79-97`. Implementaciones: `src/Mesher.cpp:47-58
(refineMesh)` y `:198-206 (generateMesh)`.

### 4.3 Estado por celda (`Quadrant`)

```cpp
enum class WindingState { Unknown, AllInside, AllOutside, Mixed };  // Quadrant.h:61

// Volume fraction (bloque nuevo)
unsigned int mSampleSize;        // Quadrant.h:189
vector<double> mWindingNumbers;   // Quadrant.h:190  ← s² doubles por celda
double mVolumeFraction;           // Quadrant.h:191  ← promedio
bool mHasVolumeFraction;          // Quadrant.h:192

// TUSQH (bloque nuevo)
WindingState mWindingState;       // Quadrant.h:195
```

### 4.4 Estado del bucle `windingSubdivide` (locales en `Mesher.cpp:1738-1793`)

```cpp
list<Quadrant> candidates, new_candidates;        // candidatos a subdividir este nivel
vector<Quadrant> processed;                       // celdas que ya son hoja TUSQH
map<unsigned int, unsigned int> idx_pos_map;      // q_id -> posición en processed
list<Point3D> new_pts;                            // puntos nuevos generados por SplitVisitor
vector<unsigned int> needs_classification;       // q_ids balanceados que requieren reclasificar
list<pair<unsigned int, unsigned int>> toBalance; // (q_id_vecino, q_id_procesado) a subdividir para balancear
unsigned int new_q_idx;                           // contador persistente de q_ids nuevos
```

---

## 5. Flujo de control (call-graph)

```
main()  src/Main.cpp:113
  │
  ├── parseo de flags (-T -N -E -M -n -s)        src/Main.cpp:193-216, 304-318, 389-392
  │
  └── mesher.generateMesh(...)  o refineMesh(...)   src/Main.cpp:441-454
        │
        │   ┌── si useTusqh (== -T):
        │   │     Mesher::windingSubdivide(...)     src/Mesher.cpp:1731-2218
        │   │       │
        │   │       │   por cada depth en [0, maxDepth):
        │   │       │     para cada candidato:
        │   │       │       WindingNumberSubdivisionVisitor::visit   src/Visitors/WindingNumberSubdivisionVisitor.cpp:90-180
        │   │       │         ├── (si -E) IntersectionsVisitor  → fuerza subdivisión si intersecta arista
        │   │       │         └── siempre: s × s × Polyline::windingNumber     src/Polyline.cpp:212-238
        │   │       │                clasifica en WindingState (AllInside/AllOutside/Mixed)
        │   │       │     las Mixed → SplitVisitor → 4 hijos
        │   │       │     balanceo one-irregular (toBalance)
        │   │       │
        │   │       ├── post-proceso: descarta AllOutside (VF=0)       src/Mesher.cpp:1966-2001
        │   │       │
        │   │       └── si tusqhExtraResolveDepth>0 (== -M):
        │   │             resolve pass (TUSQH puro, sin balance)        src/Mesher.cpp:2023-2202
        │   │
        │   └── si no useTusqh:
        │         generateQuadtreeMesh / splitQuadrants  (pipeline clásico, sin cambios)
        │
        └── Mesher::computeVolumeFractions(input, mSampleSize)    src/Mesher.cpp:2223-2247  ← SIEMPRE, sea TUSQH o clásico
              │
              │   para cada quadrant q:
              │     WindingNumberVisitor::visit(q)               src/Visitors/WindingNumberVisitor.cpp:65-90
              │       └── s × s × Polyline::windingNumber
              │             Quadrant::computeVolumeFraction(wn)  src/Quadrant.cpp:179-188
              │
              └── si useSubgrid (== -J):
                    Mesher::computeSubcellVolumeFractions(...)   src/Mesher.cpp:2876-2943
                      └── popula mEdgeSubcellVF / mVertexSubcellVF
                    Mesher::resolveArchipelagos(...)            src/Mesher.cpp:3286-3587
                      │
                      │   bucle iter (hasta convergencia):
                      │     computeSubcellVolumeFractions         (re-poblar tras cada split)
                      │     BFS componentes (MapEdges)
                      │     para cada arista info[2]==max con VF>=τ:
                      │       computeExteriorDirection            src/Mesher.cpp:2969-2990
                      │       isEdgeOnDomainBoundary               src/Mesher.cpp:2992-3050 (Issue #1)
                      │       bridgeSplitAtEdge                    src/Mesher.cpp:2969-3145 (Issue #2, #3)
                      │     drop componentes con < L celdas
                      │
                      └── guard: si Quadrants.empty() ⇒ skip saveOutputMesh  src/Mesher.cpp:3555..3568
              │
              └── dump VTK (si VTKOUT==true):
                    VolumeFractionVTKWriter::writeQuadTreeWithVF  →  volume_fraction_debug.vtk
                    VolumeFractionVTKWriter::writeVFHeatmap        →  volume_fraction_debug_samples.vtk
```

Salidas VTK adicionales cuando `useTusqh`:

- `Mesher.cpp:2209-2217` → `<name>_tusqh.vtk` y `<name>_winding_state.vtk`
  (vía `VolumeFractionVTKWriter::writeWindingState`).

Salidas VTK adicionales cuando `useSubgrid` (`-J`):

- `<name>_subcell_vertex.vtk` — vértices con `mVertexSubcellVF`
  (`Mesher.cpp:2940`).
- `<name>_subcell_edge.vtk` — aristas con `mEdgeSubcellVF`
  (`Mesher.cpp:2943`).
- `<name>_postarchipelago.vtk` — quadtree tras `resolveArchipelagos`
  (sólo si `Quadrants` no quedó vacío; ver §3.11.6).

---

## 6. Salidas esperadas

Cuando se invoca con `-T -N s -E -M d -n s` y `VTKOUT==true` (definido en
`src/Mesher.h:66`):

- `<name>_grid.vtk` — quadtree inicial.
- `<name>_tusqh.vtk` — quadtree tras TUSQH.
- `<name>_winding_state.vtk` — clasificación TUSQH por celda.
- `volume_fraction_debug.vtk` — quadtree con VF por celda.
- `volume_fraction_debug_samples.vtk` — heat map de los puntos s × s.

---

## 7. Notas para un merge entre ramas

### 7.1 Archivos tocados por este feature

Lista exhaustiva de archivos **nuevos o significativamente modificados**:

```
src/Main.cpp                                  (parseo -T/-N/-E/-M/-n)
src/Mesher.h                                  (firmas + windingSubdivide + computeVolumeFractions)
src/Mesher.cpp                                (orquestación + bucle TUSQH + resolve pass + computeVolumeFractions)
src/Polyline.h / src/Polyline.cpp             (probablemente sin cambios; primitive ya estaba)
src/Quadrant.h / src/Quadrant.cpp             (WindingState, mSampleSize, mWindingNumbers, mVolumeFraction, mWindingState, getSamplePoint, computeVolumeFraction)
src/Visitors/WindingNumberSubdivisionVisitor.h / .cpp   (NUEVO)
src/Visitors/WindingNumberVisitor.h / .cpp              (NUEVO)
src/Visualization/VolumeFractionVTKWriter.h / .cpp      (NUEVO)
```

Archivos **que no deben tocarse** al mergear:

- `src/Visitors/SplitVisitor.*` (reusado tal cual).
- `src/Visitors/IntersectionsVisitor.*` (sólo se le añade un uso, no cambia).
- `src/RefinementInputSurfaceRegion.*` (es lo que activa `-s`, independiente).
- Resto del pipeline clásico (`generateQuadtreeMesh`, `splitQuadrants`,
  `removeOnSurfaceSafe`, etc.).

### 7.2 Conflictos previsibles

| Zona                                    | Tipo de conflicto esperado                                                                                            |
|-----------------------------------------|---|
| `src/Main.cpp` (flags)                  | Si la otra rama añade flags nuevos o reordena el `switch`. Mantener la lista `inout = true` (`Main.cpp:169-185`) sincronizada. |
| `src/Mesher.h` (firmas)                 | Si la otra rama añade parámetros a `generateMesh`/`refineMesh`, mantener el orden `sampleSize → decoration → useTusqh → tusqhSampleSize → refineOnEdgeIntersect → tusqhExtraResolveDepth`. |
| `src/Mesher.cpp` (bucle TUSQH)          | El bloque `1731-2218` es grande; cualquier refactor en `SplitVisitor` o balanceo toca `1876-1927`. Resolver manualmente. |
| `src/Quadrant.h` (campos nuevos)        | Los miembros `mSampleSize/mWindingNumbers/mVolumeFraction/mHasVolumeFraction/mWindingState` deben ir al final de la sección `protected` para no romper layout. |
| `src/Visitors/WindingNumberVisitor.cpp` / `WindingNumberSubdivisionVisitor.cpp` | Probablemente no haya conflicto (archivos nuevos), pero respetar la firma de `Quadrant::getSamplePoint`. |
| `src/Visualization/VolumeFractionVTKWriter.cpp` | Idem, archivos nuevos. Mantener el contrato de los nombres de archivo de salida (ver §6). |

### 7.3 Verificación post-merge

Tras mergear, ejecutar (los binarios ya están en `build/`):

```bash
# Smoke test sin TUSQH (debe dar el mismo resultado que antes):
./build/mesher_roi -d data/<input>.d -u out_classic -a 4 -n 2

# Smoke test con TUSQH puro (debe activar windingSubdivide y resolve=0):
./build/mesher_roi -d data/<input>.d -u out_tusqh -a 4 -T -N 2 -E -n 2 -e

# Smoke test con resolve pass:
./build/mesher_roi -d data/<input>.d -u out_tusqh_m2 -a 4 -T -N 2 -E -M 2 -n 2 -e

# Verificación de regresión del primitive (Main.cpp:217-240, caso -t):
./build/mesher_roi -d data/<input>.d -t
```

Y comparar las salidas `volume_fraction_debug.vtk` y `_winding_state.vtk`
para confirmar que la clasificación TUSQH no se ha desplazado.

### 7.4 Reglas de oro al tocar este código

1. **Una sola fuente de winding number.** Si necesitas `wn(P)`, llama
   `Polyline::windingNumber(P)`. No lo recalcules en los visitantes.
2. **`WindingState` se escribe, no se infiere en sitios distintos.** Los
   únicos lugares que deben llamar a `q->setWindingState(...)` son
   `WindingNumberSubdivisionVisitor::visit` y, opcionalmente,
   `WindingNumberVisitor::visit` (no lo hace hoy).
3. **El balanceo de `windingSubdivide` usa el mismo `toBalance` que el
   pipeline clásico.** Si modificas `SplitVisitor` revisa
   `Mesher.cpp:1876-1927`.
4. **El resolve pass (`-M`) es deliberadamente sin balanceo.**
   `Mesher.cpp:2061-2070` inicializa `emptyToBalance` y `emptyIdxPosMap`
   para que `SplitVisitor` no haga balance. Conservar este contrato.
5. **`-n` y `-N` son independientes.** No compartir ni reutilizar el
   `mSampleSize` de `Quadrant` entre los dos visitantes; cada uno fija
   su propio `s` en la celda.

---

## 8. Glosario rápido

| Término                     | Definición                                                                                                  |
|-----------------------------|---|
| `WindingState`              | Enum por celda: `Unknown / AllInside / AllOutside / Mixed`. Ver `src/Quadrant.h:61`.                       |
| `wn(P)`                     | Winding number del punto `P` respecto al polígono de entrada. Ver `src/Polyline.cpp:217`.                 |
| `mSampleSize`               | Tamaño `s` de la grilla de muestras. Distinto entre `-n` y `-N`.                                            |
| `mWindingNumbers`           | Vector de `s²` doubles con los `wn` de los puntos de muestra de la celda.                                   |
| `mVolumeFraction`           | Promedio de los `mWindingNumbers` (entre 0 y 1).                                                            |
| `candidates` (TUSQH)        | Lista de celdas que aún pueden subdividirse en el nivel actual.                                             |
| `processed` (TUSQH)         | Celdas que ya son hoja en el nivel actual (AllInside/AllOutside o agotadas).                                |
| `toBalance` (TUSQH)         | Pares `(vecino, procesado)` que requieren subdivisión para mantener el quadtree one-irregular.              |
| `resolve pass`              | Iteraciones extra (controladas por `-M`) que subdividen las celdas `Mixed` que llegaron a `qrl == maxDepth`. |
| `WindingNumberVisitor`      | Visitor que **sólo** computa `mWindingNumbers`/`mVolumeFraction`.                                            |
| `WindingNumberSubdivisionVisitor` | Visitor que **decide** subdividir según el criterio TUSQH.                                            |

---

## 9. Historial relevante (para contexto de merge)

```
0940c2c  using winding numbers for quadtree
5f7af55  Implement Tusqh-style s×s grid winding numbers with configurable size
8c46d1d  Add sub-element centroids with individual winding numbers to debug output
ff74e08  Add centroids and individual winding numbers to debug output
2766c2d  Fix: Direct volume fraction calculation from sub_elements centroids
ad4603b  Fix: Implement proper winding number inheritance in QuadTree
8433d88  WIP: Add volume fraction computation with winding numbers
```

Si se hace un merge con `master` o `develop`, estos commits son los que
introdujeron el feature en esta rama.
