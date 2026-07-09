# Guía de implementación: TUSQH, Winding Numbers y Flags

> **Propósito:** Este documento es la guía de referencia para entender **dónde**
> y **cómo** está implementado en el código todo lo relativo a *winding numbers*,
> subdivisión estilo **TUSQH**, *volume fractions* y los flags `-T -N -E -M -n -s`
> del ejecutable. Está pensado como apoyo para futuros *merges* entre ramas
> (e.g. `develop-felipe` ↔ `develop` ↔ `master`).
>
> Auditado sobre la rama `develop-felipe` (HEAD `0940c2c using winding numbers
> for quadtree`). Las referencias `archivo:línea` son absolutas dentro del repo
> `MixedQuadTree`.

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
              └── dump VTK (si VTKOUT==true):
                    VolumeFractionVTKWriter::writeQuadTreeWithVF  →  volume_fraction_debug.vtk
                    VolumeFractionVTKWriter::writeVFHeatmap        →  volume_fraction_debug_samples.vtk
```

Salidas VTK adicionales cuando `useTusqh`:

- `Mesher.cpp:2209-2217` → `<name>_tusqh.vtk` y `<name>_winding_state.vtk`
  (vía `VolumeFractionVTKWriter::writeWindingState`).

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
