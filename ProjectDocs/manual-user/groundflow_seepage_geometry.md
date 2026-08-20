# groundflow_seepage_geometry

## Description

Specifies an edge of the groundflow domain for which the groundwater is only
allowed to flow outwards of the domain; flow into the domain is not allowed on
that edge. The geometrical entity should be specified such that the normal of
the geometry points outwards the material (so outwards the groundflow domain).

This option comes in handy when the point of groundwater flow exit is not known
in advance of the calculation; it will be a result of the calculation instead.

Typical usage combines the seepage edge with a zero-pressure boundary condition
(free air):

```
groundflow_seepage_geometry 10 -geometry_line 100
bounda_dof 20 -geometry_line 100 -total_pressure
bounda_time 20 0.0
```

The total pressure (pore pressure) is set to 0 on the geometry line 100, to
account for free air at that edge. Since at that edge water cannot enter the
domain the seepage option is applied to that edge. The result of these combined
options is that on nodes with outward flow a total pressure 0 boundary
condition is imposed, whereas on other nodes no boundary condition is imposed
(so that the flow is 0 at those nodes). The transition point between these
outflow nodes and nodes with zero flow will be found automatically as a result
of the calculation.

## Usage

```
groundflow_seepage_geometry <index> <geometry_item_name> <geometry_item_index>
```

## Parameters

| Parameter | Meaning                                                         |
|-----------|-----------------------------------------------------------------|
| `index`   | Record index.                                                    |
| `geometry_item_name` | Geometry type, e.g. `-geometry_line`.                  |
| `geometry_item_index` | Geometry index.                                        |

## Example

```
groundflow_seepage_geometry 10 -geometry_line 100
bounda_dof 20 -geometry_line 100 -pres
bounda_time 20 0.0
```

The bottom edge (line 100) is a seepage edge: water can only leave the domain
there, and the exit point is determined automatically.
