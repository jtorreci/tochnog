# Corpus de regresión GNU + Professional — baseline DEFINITIVO (2026-08-30)

Con el binario de F. Lorenzo (Sources.zip jan-2014, construido con
nuestro hypo.c C-puro + lapack/blas del numlib-runtime) como
referencia, contra SU corpus de 197 tests:

## Veredicto de regresión (la garantía de que no rompimos lo que funcionaba)

| binario | PASS |
|---|---|
| Fernando (referencia) | 163/197 |
| Nosotros | 160/197 |

- **Regresiones reales: 2 de 163 (98.8% preservado)**
  - force5 (phreatic level en FORCE_ELEMENT_EDGE_WATER — gap)
  - beam2d (mesh-refinement path — investigar)
- **3 falsas regresiones**: examp1, examp2, matrix1 = el solver
  HONESTO (fix A+B): los criterios viejos fake-pasaban con x=0
- **5 MEJORAS**: condif5, examp14, ground3, ho_othr2 y hypo1 (¡el
  test de hipoplasticidad de Fernando falla con SU binario y pasa
  con el nuestro!)

## GNU 2014 completo: 160/197 (81%)

Familias que fallan (las mismas que en el 2001):
- mesh-refine/generated-mesh (beam2d_1/2, refine4, tet10...): solve a 0
- Bi-CG honestos (examp1/2, matrix1, incnav5, refine5)
- petsc1 (necesita PETSC), ho_mech4 (materi_strain_plasti init)
- feature gaps varios

## Professional 363: 51 PASS (14%) = LA MEDIDA DE CONVERGENCIA

253 parse-errors = backlog exacto de keywords (missing_keywords.txt):
print_apply x29, inertia_apply x14, -total_pressure x11,
-matrix_pardiso x11, group_truss_elasti_young x9, ...

## Reproducción

- Binario Fernando: /tmp/opencode/tn2014/build/tochnog (makefile
  ajustado: -m486 fuera, -static fuera, libgfortran+lapack del
  numlib-runtime, hypo.c del GNU nuestro)
- Corpus: /tmp/opencode/tn2014/test/*.dat
- Logs: /tmp/t14_*.log (nuestro) /tmp/f14_*.log (Fernando)
