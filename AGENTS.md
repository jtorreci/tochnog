# Tochnog — Project conventions

## Documentation of new features (MANDATORY)

Every time a feature is implemented (keyword, data record, or behavior change),
you MUST document it in TWO places, in the same work unit as the feature:

1. `ProjectDocs/manual-user/<feature>.md` — USER manual.
   Content: what the feature does, why it is useful, input syntax, parameters
   with their meaning, and a minimal example input.

2. `ProjectDocs/manual-developer/<feature>.md` — DEVELOPER manual.
   Content: where it was implemented (file + function), implementation details,
   external features/libraries used, hardcoded parameters, and pending
   refactorings (cryptic enums, hardcoded data).

Also update the index files:
- `ProjectDocs/manual-user/README.md`
- `ProjectDocs/manual-developer/README.md`

Rules:
- Use the exact feature name as the .md filename (lowercase, snake_case).
- Write in English.
- Commit the documentation together with the feature implementation.
- If a feature is only partially implemented, mark it clearly and note what
  remains (user) and which parts are pending (developer).

## Convergence tracking (MANDATORY)

The project's goal is convergence with Tochnog Professional. Track every
feature in `ProjectDocs/SEGUIMIENTO-CONVERGENCIA.md`:

- When implementing a feature: update its line to `[x]` with the commit
  hash and date, and add a row to the "registro de verificación" table.
- When deciding NOT to implement a feature: add it to the "Descarte de
  features" table with the reason, and mention it in the "Differences with
  the Professional version" section of the relevant manual.
- Verification rule: before implementing a feature, search for a usage
  example — in the Professional manual, in its repo/Google Drive (links on
  the Professional website), or in the sfnet test suite. Verify against
  that example if found; otherwise build a dedicated test. Record how the
  feature was verified.

## Build notes

- Only clean builds are reliable on the current filesystem (`/mnt/z`).
  Incremental builds mix old/new enums and corrupt the binary.
- Linking to `/tmp` and copying the binary to `build/` avoids I/O failures.
- `tochnog.h` and `tochnog-mod.h` must stay in sync (same enum, same order).
