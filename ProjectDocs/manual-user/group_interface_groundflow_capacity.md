# group_interface_groundflow_capacity

## Description

Specifies the capacity for interface elements in a ground flow analysis. The
interface stores water when its pore pressure changes (storage term on the
pressure dofs of the interface nodes).

## Usage

```
group_interface_groundflow_capacity <element_group> <C>
```

## Parameters

| Parameter | Meaning                                                |
|-----------|--------------------------------------------------------|
| `element_group` | Element group, see `element_group`.            |
| `C`        | Storage capacity of the interface.                     |

## Example

```
group_type 10  -groundflow
group_interface 10  -yes
group_interface_groundflow_capacity 10  1.0
```
