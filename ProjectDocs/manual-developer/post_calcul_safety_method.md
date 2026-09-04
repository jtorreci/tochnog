# post_calcul_safety_method — developer

Enum + registration (INTEGER, 1 value, no_index 1, class POST). Pure-name
selectors: PRIVAL ("prival") pre-existed; VERTICAL/GLOBAL ("vertical"/"global")
added as pure name entries. Operator values SAFETY_PIPING/SAFETY_LIFTING
registered as INTEGER pure names so the post_calcul record parses the
-safety_piping/-safety_lifting operators. Label naming measured against the
Professional .dbs on ground15/16: vertical -> "safety_piping"/"safety_lifting"
(1 value), prival -> "safety_piping_prival_0..2" (3 values), global ->
"safety_piping_global_x/y/z" (3 values). Computation branches in calcul.cc
PENDING (formulas: lifting = (sigma_i+p_total)/p_total, piping =
(sigma_i+p_dynamic)/p_dynamic, calibrated on the corpus values 2.0/2.3333).
Blocked upstream by the groundflow pressure conventions (-topres bounda).
