# groundflow_pressure_factor

## Description

Multiplicative factor applied to the pore pressure when computing the total
stress in the `materi_stress` model:

```
sigtotal = sig_efectivo + factor * pres
```

The effective stress computed by the material model is left untouched; only the
pore-pressure contribution added to the total stress is scaled. With the
default value `1` the behaviour is identical to not defining the record.

Useful for consolidation/stability problems, for example to scale the safety
factor on the pore pressure when assessing stability.

Requires the `groundflow_pressure`, `materi_stress` and `materi_velocity`
models.

## Usage

```
groundflow_pressure_factor <factor>
```

The record is `no_index`: a single value applies to the whole model.

## Parameters

| Parameter | Meaning                                                        |
|-----------|----------------------------------------------------------------|
| `factor`  | Multiplier of the pore pressure added to the total stress. Default `1`. |

The value multiplies the pore pressure `pres` before it is added to the total
stress:

```
total_new_sig[diagonal] += factor * pres
```

## Example

```
groundflow_pressure_factor 2.
```

The pore pressure contribution to the total stress is doubled, e.g. to run a
stability analysis with a safety factor of 2 on the pore pressure.
