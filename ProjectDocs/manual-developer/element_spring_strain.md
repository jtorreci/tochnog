# element_spring_strain — developer

Enum + registration: DOUBLE, data_length 1, version_all 1 (written per step at
VERSION_NEW like element_spring_force), print_only 1, class ELEMENT,
data_required ELEMENT. Pre-allocated in top.cc inside the any_spring block
(step_start). Written in spring.cc only when the nonlinear diagram is active
(currently); targets read VERSION_NORMAL after the step promotion.
