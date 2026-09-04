# post_calcul_safety_method

Determines how the hydraulic safety factors of `post_calcul -materi_stress
-safety_piping` / `-safety_lifting` are computed (manual Professional 6.919):
`-vertical` (default, the vertical normal stress; one value), `-prival` (the
three principal stresses; three values `..._prival_0..2`) or `-global` (the
three global normal stresses; three values `..._global_x/y/z`).

Formulas (measured against the Professional .dbs on ground15/16 of the corpus):
`safety_lifting = (sigma_i + p_total)/p_total` and
`safety_piping = (sigma_i + p_dynamic)/p_dynamic`, with tochnog sign
conventions (compression negative); p_total/p_dynamic follow the groundflow
total/static/dynamic pressure split.

REGISTERED in this batch (values and label naming verified against the
Professional: `safety_piping_prival_0`, `safety_lifting_global_y`, ...). The
computation branches in calcul.cc are PENDING together with the groundflow
pressure conventions (bounda_dof -topres) that block the family.
