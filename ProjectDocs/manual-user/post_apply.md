# post_apply

## Description

Global switch of the post-processing commands (manual Professional
6.900): `post_apply -no` prevents the post-processing commands
(records with "post" in the name) from being evaluated; only the
`post_node_rhside_ratio` is evaluated always. Default is `-yes`.

## Usage

```
post_apply [-yes | -no]
```

## Status

Registered (INTEGER, no_index) and parsed. The `-no` gating of the
per-step post evaluation is PENDING (no corpus test uses `-no`; the
distri3 test only sets `-yes`, the default). `print_apply` (the print
counterpart, 6.968) was already implemented.
