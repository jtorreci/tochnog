# Corpus de regresión GNU + Professional — baseline (2026-08-30)

Tras el incidente del rm (ver SEGUIMIENTO), el corpus de tests se
reconstruye desde las DOS distribuciones originales:

## GNU 2001 (202 tests, /tmp/opencode/tn2001 + SourceForge)

- **164 PASS (81%)** con nuestro binario actual — el motor está sano
- 38 FAIL:
  - ~16 camino de refinado/malla-generada (beam2d familia): solve a 0
    (machinery de mesh-change, NO tocado por nosotros; investigación
    aparte)
  - 6 Bi-CG breakdown: CONSECUENCIA ESPERADA del solver honesto (fix
    A+B — los criterios viejos "pasaban" con x=0 basura; documentado
    en DIAG-SOLVE-MIXTO): examp1 examp2 matrix1 incnav5 refine5
  - petsc1: necesita PETSC (compilación)
  - resto: gaps de features / sintaxis 2001

## Professional 363 tests vs binario GNU — LA MEDIDA DE CONVERGENCIA

- **51 PASS (14%)**
- **253 PARSE-ERROR (70%)** = el backlog exacto de convergencia
- 59 RUN-FAIL (16%)

### Top keywords faltantes (de los parse errors)

29 print_apply · 14 inertia_apply · 11 -total_pressure ·
11 -matrix_pardiso · 9 group_truss_elasti_young ·
7 groundflow_phreatic_level · 5 -quad6 · 5 -contact_spring2 ·
4 -truss_beam · 4 -size_dev · 4 geometry_factor ·
3 post_element_force · 2 -updated_area · 2 slide_plasti_friction ·
2 -quad8 · 2 processors · 2 -nory_sig ·
2 group_materi_undrained_capacity · 2 geometry_element_group ·
2 -contact_spring (+ ~30 más en missing_keywords.txt)

## Uso como suite de regresión

- GNU: `for f in /tmp/opencode/tn2001/tochnog/test/*.dat; do ...`
  (script en la memoria de sesión)
- Professional: idem sobre test/{validation,tutorial,other}
- El 14%→100% del corpus Professional ES la convergencia
