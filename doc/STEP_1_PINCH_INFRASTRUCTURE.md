# Step 1 — Pinch detection infrastructure

**Status:** Pendiente
**Date:** 2026-07-16
**Branch:** `feature/paper-faithful-bridge`
**Depends on:** Step 0 (commit baseline) ✅

---

## Estado actual del proyecto (post-Fase 0)

### Git
- HEAD: `0912da9 Paper-faithful TUSQH bridge-joining + sub-cell VF (steps 1-7)`
- Working tree: limpio
- Rama: `feature/paper-faithful-bridge` basada en `develop-felipe@d174e6c`

### Verificación reciente
- `cmake --build build` → clean (100% Built target mesher_roi)
- `python3 scripts/test_bridge_clean_info.py --mesher build/mesher_roi` → **6/6 pass**

### Funcionalidad del proyecto (paper §3.1-§3.4e)
- ✅ §3.1 Volume Fractions (grilla + sampling)
- ✅ §3.2 Persistence Parameter (umbral)
- ✅ §3.3 Sub-cell VFs (vértices + aristas)
- ✅ §3.4e Archipelago resolution (paper-faithful)
- ❌ §3.4a Pinch detection (vértices, 2×2)
- ❌ §3.4c Templates separating (Fig. 11)
- ❌ §3.4d Consistencia entre pinches adyacentes (Fig. 13)

---

## Objetivo de Fase 1

Crear las **clases base** y **estructuras de datos** que las Fases 2-5 usarán para implementar la detección y resolución de pinches del paper §3.4a-d.

### Principios de diseño

1. **Cero impacto en el pipeline actual.** El flag `-P` será default 0 (off). Sin `-P`, el comportamiento debe ser idéntico al actual.
2. **API mínima.** Cada clase expone sólo lo necesario para Fases 2-6.
3. **Estructuras reutilizables.** `PinchCase` enum, `VertexPinchInfo` struct, etc.
4. **Documentación inline.** Cada nueva función tiene un Doxygen corto.

### Archivos a crear

| Archivo | Líneas (aprox) | Propósito |
|---------|----------------|-----------|
| `src/PinchCase.h` | 60 | `enum class PinchCase` + struct `VertexPinchInfo` + struct `PinchConfig` |
| `src/Visitors/PinchVisitor.h` | 90 | Clase abstracta para visitar vértices (estilo Visitor existente) |
| `src/Visitors/PinchVisitor.cpp` | 50 | Constructor trivial + setters |
| `src/PinchDetector.h` | 100 | API de detección: `detectAtVertex`, `detectAll` |
| `src/PinchDetector.cpp` | 80 | Stubs que serán implementados en Fase 2 |

### Archivos a modificar

| Archivo | Líneas añadidas | Cambios |
|---------|----------------|---------|
| `src/Main.cpp` | ~20 | Variable `bool usePinchDetection = false`; parseo de `-P` (acepta 0/1/2) |
| `src/MeshPoint.h` | ~15 | Campos `mPinchCase`, `mIsPinchVertex` (sin uso aún) |
| `src/MeshPoint.cpp` | ~5 | Inicialización de los nuevos campos |
| `src/CMakeLists.txt` | ~5 | Añadir `PinchDetector` a SOURCE_FILES |

### Sin tocar
- `src/Mesher.h` / `src/Mesher.cpp` — los pinch handlers se añadirán en Fases 3, 4, 5.
- `src/Visitors/SplitVisitor.cpp` — la 1-a-5 ya existe como `bridgeSplitAtEdge`.
- `src/resolveArchipelagos` — la integración es Fase 6.

---

## Diseño de las estructuras

### `PinchCase.h`

```cpp
namespace Clobscode {

// Paper §3.4 / Fig. 9: classification of a grid vertex regarding
// its neighborhood in the cubical complex.
enum class PinchCase {
    None,                    // no pinch detected
    Vertex_2x2_Chess,        // paper Fig. 9a: 2 diagonal interior + 2 diagonal exterior (regular 2x2 case)
    Vertex_Hanging_2,        // paper Fig. 9c/d: only 2 incident quads (T-junction)
    Vertex_Hanging_3,        // paper Fig. 9e: 3 incident quads (T-junction corner)
    Vertex_Corner_1          // paper Fig. 9l: 1 incident quad (domain corner)
};

// Resolution mode required by the paper.
enum class PinchResolution {
    None,        // not a pinch
    Connect,     // paper Fig. 12 (growing template, 1-to-5 keep bridge quad)
    Separate     // paper Fig. 11 (shrinking template, 1-to-5 remove child)
};

// Per-vertex pinch state.
struct VertexPinchInfo {
    PinchCase pinCase = PinchCase::None;
    PinchResolution resolution = PinchResolution::None;
    // Indices into the Quadrants vector of the incident quads.
    // Up to 4 quads for a regular 2x2 vertex.
    std::vector<unsigned int> incidentQuads;
    // For Vertex_2x2_Chess: indices of the two interior diagonal quads.
    unsigned int interiorDiagA = 0;
    unsigned int interiorDiagB = 0;
    // Sub-cell VF at the vertex (for resolution decision).
    double subcellVF = 0.0;
};

} // namespace Clobscode
```

### `PinchVisitor`

Patrón Visitor existente ya tiene una clase base (`src/Visitors/Visitor.h`). Como pinch opera sobre **vértices** (no quads), no encaja en el patrón existente. Se hace un visitor standalone:

```cpp
namespace Clobscode {

class PinchVisitor {
public:
    PinchVisitor();
    virtual ~PinchVisitor();

    void setPoints(const std::vector<MeshPoint> *mp);
    void setQuadrants(const std::vector<Quadrant> *q);
    void setMapEdges(const std::map<QuadEdge, EdgeInfo> *me);
    void setJoinThreshold(double t);

    // Called for each vertex of the quadtree.
    // Returns true if the visitor wants to continue iterating.
    virtual bool visit(unsigned int vertexId,
                       const std::vector<unsigned int> &incidentQuads) = 0;

protected:
    const std::vector<MeshPoint> *mPoints = nullptr;
    const std::vector<Quadrant> *mQuadrants = nullptr;
    const std::map<QuadEdge, EdgeInfo> *mMapEdges = nullptr;
    double mJoinThreshold = 0.5;
};

} // namespace Clobscode
```

### `PinchDetector`

Stubs en Fase 1; implementación en Fase 2.

```cpp
namespace Clobscode {

class PinchDetector {
public:
    PinchDetector();

    // Detects pinch cases in the quadtree.
    // Output: mPinchInfo[i] = pinch state for vertex i.
    // Returns the number of pinches detected.
    unsigned int detectAll();

    // Read-only access to per-vertex pinch state.
    const std::vector<VertexPinchInfo> &getPinchInfo() const { return mPinchInfo; }
    unsigned int getPinchCount() const;

private:
    // Phase 2 will implement this.
    PinchCase detectAtVertex(unsigned int vertexId,
                              const std::vector<unsigned int> &incidentQuads);

    std::vector<VertexPinchInfo> mPinchInfo;
};

} // namespace Clobscode
```

---

## Comandos de Fase 1

### Pre-modificación

```bash
# Confirmar que estamos en la rama correcta
git branch --show-current
# Esperado: feature/paper-faithful-bridge

# Confirmar que build está limpio
cmake --build build
```

### Creación de archivos

1. `src/PinchCase.h` — enum + structs
2. `src/Visitors/PinchVisitor.h` — header del visitor
3. `src/Visitors/PinchVisitor.cpp` — implementación trivial
4. `src/PinchDetector.h` — header del detector
5. `src/PinchDetector.cpp` — stubs

### Modificación de archivos existentes

1. `src/Main.cpp` — añadir variables `usePinchDetection`, parseo `-P`
2. `src/MeshPoint.h` — añadir campos pinch
3. `src/MeshPoint.cpp` — inicializar campos pinch
4. `src/CMakeLists.txt` — añadir PinchDetector a SOURCE_FILES

### Verificación post-Fase 1

```bash
cmake --build build
# Esperado: clean compile, no warnings nuevos

python3 scripts/test_bridge_clean_info.py --mesher build/mesher_roi
# Esperado: 6/6 pass (comportamiento idéntico, -P default off)

./build/mesher_roi -p data/unit_square.poly -u out -a 2 -P 1 -e
# Esperado: funciona sin errores (aunque todavía no hace nada con -P)
```

---

## Riesgos de Fase 1

| Riesgo | Mitigación |
|--------|-----------|
| Compilación rota por signatures nuevas en MeshPoint | Los campos nuevos son triviales (PinchCase default None, bool default false) |
| `-P` flag entra en conflicto con flags existentes | Verificar primero que `-P` no esté usado. (-P no está en `Main.cpp` actualmente.) |
| PinchVisitor no encaja con patrón Visitor existente | Se documenta como "visitor-style standalone" (no extiende `Visitor`) |

---

## Verificación esperada

- **Antes:** Comportamiento idéntico (sin `-P`).
- **Después:** Sin `-P` mismo comportamiento. Con `-P 1` o `-P 2` el binario acepta el flag pero no hace nada todavía (la lógica llega en Fases 2-5).
- **Tests:** 6/6 siguen pasando.

---

## Siguiente fase

**Fase 2 — Detección de pinches 2×2** (paper §3.4 + Fig. 9a).

Implementar `PinchDetector::detectAtVertex` para reconocer el patrón "ajedrez" del paper (2 quads diagonales interior + 2 quads diagonales exterior). Sin esto, el resto del flujo no tiene pinches que procesar.
