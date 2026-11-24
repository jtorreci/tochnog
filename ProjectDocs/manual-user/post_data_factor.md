# post_data_factor

## Description

`post_data_factor index factor_0 factor_1 ...` (manual Professional 6.925): one factor per entry of `post_data` with the same index. A shorter record is an error ("post_data_factor needs a factor for each data item of post_data.").

## Usage

```
post_data_factor 0 1. 2. 10.
```

## Example

`validation-suite/test-2026/post_count_data1.dat`: three unit quad4 elements, two element groups, the right edge pulled. `post_data_factor_error1.dat`: two factors for three entries.

## Differences with the Professional version

None measured.
