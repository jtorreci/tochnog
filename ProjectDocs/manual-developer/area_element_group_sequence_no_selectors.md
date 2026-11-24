# area_element_group_sequence (without selectors)

## Implementation

- **File**: group.cc, `area_element_group_sequence()`.
- The old code errored out when neither `use_geometry` nor `use_element`
  was set. Now that case is handled:
  - `area_element_group_sequence_element[0] = -ALL` (all elements are
    candidates).
  - In the per-element selection block, when
    `!use_geometry && !use_element`, the element must currently have the
    PREVIOUS group in the sequence: the current `ELEMENT_GROUP` is
    compared against the group of the last time window strictly before
    `time_total` (the "previous group" lookup uses the same time-window
    logic as the forward mapping, with a fallback to the last group
    before the current window). If it does not match, the element is
    skipped (`ok = 0`).
- This mirrors the Professional semantics: at `time_i` the elements of
  group `group_(i-1)` become `group_i` (the previous group selects the
  elements).

## Verification

- dam_building: parses and runs the layered construction with the
  no-selector mode (the model needs ~10 minutes, above the corpus
  budget). The group migration is exercised during the layer
  construction.

## Pending

- A dedicated small test (2 groups, sequence of 2 time windows, no
  selectors, checking the element group at each time) is pending - the
  only corpus consumer (dam_building) is too slow for the corpus budget.
