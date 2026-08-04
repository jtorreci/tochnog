# include

## Descripción

`include filename` reads another data file from within the data part and
processes it exactly as if its content were typed inline. This is useful to
keep big or generated input out of the main file: meshes produced by
pre-processors, material libraries, and repeated parts of the input.

Restrictions:

- The included file must end with the `end_data` keyword.
- The included file may NOT contain another `include`.
- The included file should not contain comments.

## Uso

Place it in the data part, as a keyword line with the file name:

```
include <filename>
```

The file name is relative to the working directory. The content of the
included file is parsed by the same data-part loop, so every keyword it
contains must be a valid Tochnog keyword.

## Parámetros

| Parameter | Meaning                                                    |
|-----------|------------------------------------------------------------|
| `filename`| Name of the file to include. It is read from the working directory and its records are processed as if they were part of the main input. |

## Ejemplo

Minimal input that reads a mesh and a material library from separate files:

```
include mesh.dat
include materials.dat

control_geometry
   cartesian
control_time
   end 1.0
   dt 0.1
control_print
   history 0
   step 10
end_data
```

Where `mesh.dat` contains the `node` and `element` records and
`materials.dat` the `materi_*` records, each one closed with its own
`end_data`.
