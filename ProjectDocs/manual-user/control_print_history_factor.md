# control_print_history_factor

## Descripción

`control_print_history_factor` defines multiplication factors for the
values written by `control_print_history`. One factor is given per printed
data item. It is useful to scale output units without changing the model
input, for example converting Pa to kPa.

## Uso

Place it in the data part, following the same index convention as
`control_print_history`:

```
control_print_history_factor 0 f0 f1 ...
```

## Parámetros

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `0`       | Index of the `control_print_history` record the factors belong to. |
| `f0 f1 ...` | One multiplication factor per data item in the matching `control_print_history` record. Default is 1.0 for items without a factor. |

## Ejemplo

Print stresses in kPa (input model is in Pa):

```
control_print
   history 0
   step 10
control_print_history
   0 2 1 2 4 2 5
control_print_history_factor
   0 1.e3 1.e3 1.e3
```

The three history data items (e.g. two stress components and a strain)
are written multiplied by 1.e3.
