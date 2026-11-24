# materi_plasti_maximum_iterations

## Description

Plasticity iterates at integration point level until the stress is back on
the yield surface. `materi_plasti_maximum_iterations` (manual Professional
6.804) caps the number of those iterations for all element groups.
`group_materi_plasti_maximum_iterations` sets the cap for one element group
and wins over the global record for that group.

Raise the cap when results of a calculation with very heavy plasticity look
inaccurate. Lower it to stop slow or non-converging plasticity loops early:
when the cap is reached, the last iterate is used without an error.

## Syntax

```
materi_plasti_maximum_iterations  maximum_iterations
group_materi_plasti_maximum_iterations  index  maximum_iterations
```

| Parameter | Meaning |
|---|---|
| `maximum_iterations` | Maximum number of plastic iterations per integration point; at least 1. |
| `index` | Element group (`group_materi_plasti_maximum_iterations` only). |

Without either record the cap is 20000 (manual default), 10 for
viscoplastic groups and 2 for nonlocal viscoplastic groups. A loop that
converges earlier never reaches the cap.

## Example

Mohr-Coulomb square in compression, capped at one iteration: the return
to the yield surface is not reached and sigyy stays elastic.

```
group_materi_plasti_mohr_coul 0 0.5 1. 0.
materi_plasti_maximum_iterations 1
```

Tests: `validation-suite/test-2026/materi_plasti_maximum_iterations1`
(cap 1, sigyy -21.978), `materi_plasti_maximum_iterations2` (group 1000
wins over global 1, sigyy -3.3935) and
`materi_plasti_maximum_iterations_error1` (0 is an error);
`materi_plasti_maximum_iterations_default1` (default cap 20000 reached
on the sfnet examp16 foundation).

## Differences with the Professional version

- Default 20000 as the manual since 2026-11-06 (owner decision, BACKLOG
  C112; before: 1000).
- The Professional accepts 0 and runs as with 1; the GNU stops with an
  error (out-of-domain value).
- With the cap at 1 the Professional gives sigyy -21.9773 against the GNU
  -21.978: the iteration schemes differ in detail.
- `group_materi_plasti_maximum_iterations` is a GNU record; the
  Professional has only the global one.
