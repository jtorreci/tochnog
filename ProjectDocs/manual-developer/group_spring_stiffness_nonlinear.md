# group_spring_stiffness_nonlinear — developer

Enum + registration (DOUBLE, variable DATA_ITEM_SIZE, class SPRING,
data_required GROUP_TYPE). Consumed in spring.cc (spring()): the element
activates when EITHER the linear stiffness or the nonlinear diagram exists;
the diagram is read via db_dbl/db_len as (eps,k) pairs. The stiffness of the
increment is evaluated at the MIDPOINT total strain of the increment
(spring_strain_old + spring_strain_total)/2 so a piecewise linear diagram is
integrated EXACTLY along a linear strain path (spring6 checks F=0.5 with
tol 1e-8). Strain definition: total elongation = new_length - initial_length.
Hardcoded guards: nl_length>=4 and even; flat extrapolation outside the
diagram. ELEMENT_SPRING_STRAIN is allocated in top.cc step_start next to
ELEMENT_SPRING_FORCE (parallel loop) and PUT at VERSION_NEW per spring step.
