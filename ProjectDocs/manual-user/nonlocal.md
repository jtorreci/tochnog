# nonlocal / nonlocal_name

## Description

`nonlocal` (manual Professional 6.897) and `nonlocal_name` (6.898)
configure the nonlocal plasticity:

```
nonlocal 0.8
nonlocal_name -group_materi_plasti_mohr_coul_hardening_softening
```

- `nonlocal <radius>` is the short alias of `options_nonlocal`: the
  radius of the gaussian averaging region for the nonlocal yield rule
  `f_n` (initialized as dof by `materi_plasti_f_nonlocal`).
- `nonlocal_name <name>` selects the plasticity model that is treated
  nonlocal.

## Usage

```
materi_plasti_f
materi_plasti_f_nonlocal
end_initia

nonlocal 0.8
nonlocal_name -group_materi_plasti_mohr_coul_hardening_softening
```

## Notes

- The GNU applies the nonlocal yield-rule contribution to every
  plasticity model (plasti.cc: `f += options_nonlocal * grad(f_n)`); the
  `nonlocal_name` record is accepted but has no per-model gate
  (partial).
- slope_nonlocal_refine uses these records and now parses; the run
  hits a pre-existing buffer bug of the legacy nonlocal machinery
  (nonloc.cc, `NONLOCAL_ITEM_SIZE` overflow for dense meshes) - see
  SEGUIMIENTO.
