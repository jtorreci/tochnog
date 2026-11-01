# control_post — developer

Record `CONTROL_POST` (INTEGER, CONTROL), accepted; no consumer: the
`step_close( NO )` of a control index without timesteps already calls
`post()`, and `top()` recomputes a stale `post_element_force` there (see
[post_apply](post_apply.md)). Measured on the Professional 02-08-2026:
`control_post 2 -yes`, `-no` and absent give the same post values.
