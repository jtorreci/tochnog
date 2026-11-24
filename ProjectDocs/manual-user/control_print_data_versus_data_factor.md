# control_print_data_versus_data_factor

## Descripción

`control_print_data_versus_data_factor` defines multiplication factors for
the values written by `control_print_data_versus_data`. One factor is given
per data item. It is useful to scale output units without changing the
model input, for example converting Pa to kPa.

## Uso

Place it in the data part, following the same index convention as
`control_print_data_versus_data`:

```
control_print_data_versus_data_factor 0 f0 f1 ...
```

## Parámetros

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `0`       | Index of the `control_print_data_versus_data` record the factors belong to. |
| `f0 f1 ...` | One multiplication factor per data item in the matching `control_print_data_versus_data` record. Default is 1.0 for items without a factor. |

## Ejemplo

Print two data items scaled to kPa:

```
control_print
   history 0
   step 10
control_print_data_versus_data
   0 2 1 2 4 2 5
control_print_data_versus_data_factor
   0 1.e3 1.e3
```
