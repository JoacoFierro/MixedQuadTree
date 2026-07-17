# Step 0 — Commit baseline (paper-faithful-bridge steps 1-7)

**Status:** Pendiente
**Date:** 2026-07-16
**Branch:** `feature/paper-faithful-bridge`

---

## Estado actual del proyecto (antes del commit)

### Rama
- `feature/paper-faithful-bridge` basada en `develop-felipe@d174e6c`
- Working tree sucio: cambios sin commitear acumulados en sesiones previas
- HEAD sigue siendo `d174e6c` (no hay commits nuevos en esta rama)

### Cambios sin commitear que se van a preservar

| Archivo | Tipo | Líneas (aprox.) | Contenido |
|---------|------|-----------------|-----------|
| `src/Mesher.h` | modificado | +152 | Declaraciones nuevas (`preRefineForTusqh`, `computeSubcellVolumeFractions`, `bridgeSplitAtEdge`, `computeExteriorDirection`, `isEdgeOnDomainBoundary`, `isInteriorCell`, `resolveArchipelagos` reescrita) |
| `src/Mesher.cpp` | modificado | +1023 | Implementaciones + reescritura de `windingSubdivide` (Step 1 preservación AllOutside) + bucle bridge-joining paper-faithful (Step 3) + output filter (Step 4) + heurística pre-refine |
| `src/Main.cpp` | modificado | +57 | Flags `-T -N -E -M -J -K -F -L` |
| `src/CMakeLists.txt` | modificado | +10 | `SubgridSampler` añadido a SOURCE_FILES |
| `src/MeshPoint.h/.cpp` | modificado | +100 | Campos sub-cell VF (`mSubcellVolumeFraction`, etc.) |
| `src/Visualization/VolumeFractionVTKWriter.{h,cpp}` | modificado | +211 | Writers para `_subcell_*.vtk` |
| `src/SubgridSampler.{h,cpp}` | NUEVO | ~350 | Lógica de muestreo sub-celda |
| `src/SubcellVFData.h` | NUEVO | ~50 | POD `EdgeSubcellVFData` |
| `doc/TUSQH_WINDING_GUIDE.md` | modificado | +305 | §3.11.9 reescrito con paper-faithful |
| `doc/BUGS_FOUND.md` | NUEVO | 455 | 7 issues documentados |
| `doc/WORK_SUMMARY.md` | NUEVO | 181 | Resumen maestro + métricas finales |
| `doc/FUTURE_WORK.md` | NUEVO | 262 | TODOs P1-P5 |
| `doc/CURRENT_STATE.md` | NUEVO | 155 | Estado del branch |
| `doc/README.md` | NUEVO | 109 | Índice de docs |
| `doc/SESSION_TUSQH_BRIDGE_2026-07-15.md` | NUEVO | 320 | Log de sesión |
| `doc/STEP_0..7_*.md` | NUEVOS (8) | ~50 c/u | Reportes por paso |
| `scripts/test_bridge_clean_info.py` | NUEVO | ~250 | 6 tests de regresión |
| `scripts/analyze_output.py` | NUEVO | ~200 | Validación cuantitativa |
| `data/Agua.poly` | NUEVO | 46432 vértices | Chesapeake Bay |
| `data/tusqh_*.poly` (5 archivos) | NUEVOS | varios | Datasets de prueba |
| `build/`, `*.vtk` en raíz | NO commitear | — | Artefactos de compilación y debug |

### Funcionalidad ya implementada (paper §3.1-§3.3 + §3.4e)

✅ **§3.1 Volume Fractions:** `generateGridMesh` + `windingSubdivide` + `computeVolumeFractions`.
✅ **§3.2 Persistence Parameter:** Umbral `-F` decide retención por VF.
✅ **§3.3 Sub-cell VFs:** `computeSubcellVolumeFractions` para vértices (`mVertexSubcellVF`) y aristas (`mEdgeSubcellVF`).
✅ **§3.4e Archipelago resolution (paper-faithful):** `resolveArchipelagos` reescrita en Step 3 (BFS sobre cubical complex completo + filtros Issue #1/#2/#3).
✅ Bugs resueltos: Issues #1-#6 (Issue #7 sólo documentado).
✅ Tests: 6/6 pasando.

### Métricas finales (Chesapeake Bay, `data/Agua.poly -a 3 -T -J -K 2 -F 0.5 -L 1`)

| Métrica | Valor |
|---------|-------|
| Interior cells en output | 549 |
| Componentes en output | 110 |
| Bridges añadidos | 125 (2 iteraciones) |
| Wall-clock | ~38 s |
| AllOutside cells descartadas | 970 |

---

## Qué falta para imitar el paper (paper §3.4 completo en 2D)

### Piezas del paper NO implementadas

| § del paper | Componente | Estado |
|-------------|-----------|--------|
| 3.4a | Pinch detection (vértices, 2×2) | ❌ NO existe |
| 3.4b | Templates "connecting" (growing, Fig. 12) | ✅ vía `bridgeSplitAtEdge` |
| 3.4c | Templates "separating" (shrinking, Fig. 11) | ❌ NO existe |
| 3.4d | Consistencia entre pinches adyacentes (Fig. 13) | ❌ NO existe |

### Limitaciones residuales reconocidas (de `FUTURE_WORK.md`)

- **P2 #2.3** — 3D pinch templates (no aplica, proyecto es 2D).
- **P2 #2.4** — Optimización de BFS.
- **P2 #2.5** — Diagramas de persistencia (requiere Gudhi, fuera del alcance).
- **P1 #1.1** — Verificación visual en ParaView (responsabilidad del usuario).

---

## Plan inmediato: 9 fases

```
Fase 0  — Commit baseline (ESTE PASO)
Fase 1  — Infraestructura: PinchVisitor, PinchDetector, PinchCase, flag -P
Fase 2  — Detección de pinches 2×2 (paper §3.4 / Fig. 9a)
Fase 3  — Template "separating" / shrinking (paper Fig. 11a)
Fase 4  — Consistencia entre pinches adyacentes (paper Fig. 13)
Fase 5  — Manejo de T-junctions en pinch detection
Fase 6  — Integración con resolveArchipelagos
Fase 7  — Validación contra Chesapeake Bay (paper Figs 6, 7, 22)
Fase 8  — Documentación final (STEP_8..12, PINCH_WINDING_GUIDE)
```

### Comandos de Fase 0

```bash
cd /home/felipemarchant/Proyectos/MixedQuadTree

# 1. Verificar estado actual
git status
git diff --stat

# 2. Stagear todo lo relevante (NO build/, NO *.vtk)
git add doc/ scripts/ src/ data/

# 3. Verificar que el stage es correcto
git status
git diff --cached --stat

# 4. Commit con mensaje descriptivo
git commit -m "Paper-faithful TUSQH bridge-joining + sub-cell VF (steps 1-7)

- Implement paper §3.1-§3.3 + §3.4e (archipelago resolution)
- 6 bugs fixed (Issues #1-#6 in BUGS_FOUND.md)
- 6/6 regression tests pass
- Chesapeake Bay: 549 interior cells, 110 components, 125 bridges
- See doc/WORK_SUMMARY.md and doc/STEP_0..7_*.md"

# 5. Verificar que el working tree queda limpio
git status
```

### Verificación post-commit

```bash
# Build limpio
cmake --build build

# Tests de regresión
python3 scripts/test_bridge_clean_info.py build/mesher_roi

# Salida esperada: 6/6 pass
```

### Resultado esperado de Fase 0

- HEAD de `feature/paper-faithful-bridge` avanza a un nuevo commit (ej. `a1b2c3d`).
- Working tree queda limpio (sin cambios sin commitear).
- Binario en `build/mesher_roi` sigue funcional.
- 6/6 tests siguen pasando.

---

## Riesgos de Fase 0

| Riesgo | Mitigación |
|--------|-----------|
| Olvidar excluir `build/` o `*.vtk` del commit | Verificar con `git status` antes de commitear |
| Mensaje de commit poco claro | Usar plantilla estructurada arriba |
| Cambios sin commitear adicionales no contemplados | `git diff --stat` antes y después |

---

## Siguiente fase

**Fase 1 — Infraestructura de pinches.**

Crear las clases base `PinchVisitor`, `PinchDetector`, `PinchCase` y el flag `-P` en `Main.cpp`. NO cambia comportamiento del pipeline actual — sólo añade cimientos que las Fases 2-5 usarán.
