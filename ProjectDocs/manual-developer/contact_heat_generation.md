# contact_heat_generation (developer)

## Files and functions

- `tochnog.h` / `tochnog-mod.h` — enum `CONTACT_HEAT_GENERATION` (right
  after the legacy `CONTACT_HEATGENERATION`; headers in sync).
- `database.cc` — keyword registration after the legacy one:
  `DOUBLE_PRECISION`, `data_length = 1`, `no_index = 1`,
  `data_class = CONTACT`.
- `contact.cc` — `parallel_contact()` reads BOTH names at the single
  consumption point (where `CONTACT_HEATGENERATION` was already read):
  the buffer is zeroed first (GET_IF_EXISTS does not write missing
  records — see the ddum3 rule), then the legacy name, then the
  Professional name, so the latter wins if both exist.
- The physics is inherited from the GNU 2014 sources: `friction_energy
  = eta * friction_force * slip_size` in the SLIP branch, injected into
  `node_rhside[temp_indx/nder]` of the contacter (full amount with an
  analytic `contact_geometry` target; half to each side with element
  targets, weighted for the target nodes).

## Design note: why a dual read instead of an alias PUT

The first implementation copied the Professional record into the legacy
one (PUT) at the top of `parallel_contact()`, mirroring the
`contact_target_geometry` alias. That PUT allocates database data inside
a parallel loop and aborts with "Data is allocated in a parallel loop"
(`parallel_contact` runs through `parallel_sys_routine`). GET_IF_EXISTS
reads never allocate, so the dual read at the consumption point is safe.
(The `contact_target_geometry` alias PUTs share this latent constraint;
they are only executed when the base records are absent.)

## Verification

- Test `contact_heatgen` (suite 62/62): model with `condif_temperature`
  + friction + slip under the Professional keyword name; runs stable,
  mechanics in the blocked regime (disy target).
- Value flow (documented, outside the suite): with `TOCHNOG_DEBUG=yes`
  the internal `friction_energy` scales exactly with the factor — 0.5 ->
  1.31194, 1.0 -> 2.62387 in the reference probe — proving the keyword
  value drives the physics.
- Known limitation (experimental algorithm): the accumulated node
  temperature stays ~0 in the reference model because the final
  equilibrium iteration returns to STICK (node_rhside is rebuilt each
  iteration; only the final iteration enters the solve). End-to-end
  temperature rise needs a persistently slipping contact.
