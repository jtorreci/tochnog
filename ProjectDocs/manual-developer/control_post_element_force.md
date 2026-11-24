# control_post_element_force — developer

Record `CONTROL_POST_ELEMENT_FORCE` (INTEGER, CONTROL). Consumed in the
timestep path of `top()` (`top.cc`) around `post_element_force_calculate()`
with `control_materi_gate_off()` (control only, no global form);
`post_element_force_stale` set when skipped, see
[post_apply](post_apply.md).
