# volume_factor_x — developer notes

Sprint 12 lot 3 (manual Professional 6.1094). Enum `VOLUME_FACTOR_X`,
DOUBLE, variable length (`data_length = DATA_ITEM_SIZE`,
`fixed_length = 0`), `no_index = 1`, `data_class = VOLUME` (same class
as `volume_factor`).

## Implementation

In `volume.cc::volume_factor()`, after the polynomial factor and
before the group factor (all three MULTIPLY):

```c
if ( db_active_index( VOLUME_FACTOR_X, 0, VERSION_NORMAL ) ) {
  length = db_len( VOLUME_FACTOR_X, 0, VERSION_NORMAL );
  if ( length>=3 && length%2==1 ) {
    x = coord[0];
    factor = 1.;
    if ( x>=ptr[0] && x<ptr[length-1] ) {   // inside [x0, xn)
      for ( j=0; j<length-1; j+=2 )
        if ( x>=ptr[j] ) factor = ptr[j+1];
        else break;
    }
    volfac *= factor;
  }
  else warn (odd number of values required, record neglected)
}
```

Keyed on the element INTEGRATION POINT x (`coord_ip` - the same
coordinate the polynomial factor uses), so it works for any ndim.
The right-of-last-x region is factor 1 (manual: "right from xn the
factor is 1 again") - a plain last-interval-wins loop WITHOUT the
x<xn guard would give the last factor there (bug avoided).

## Verification

`tslv_vfx`: 8x1 quad4 chain, E=1000, N=2e-2, `volume_factor_x 4. 2. 8.`
- sigxx 2e-2 / 1e-2 per half (targets), tip velx 1.15e-3 vs the
analytic 1.2e-3 (Poisson band), A/B without the record 1.59e-3 (the
1.6e-3 analytic).

## Gotchas

- An EVEN number of values warns and neglects the record (the x/factor
alternation requires odd).
