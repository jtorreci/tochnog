# geometry_element_group — developer

Enum + registration (INTEGER, variable DATA_ITEM_SIZE, class GEOMETRY).
Consumed inside geometry() (geometry.cc) at the end of the node tests: when a
GEOMETRY_ELEMENT_GROUP record exists for the geometry index being tested and
the node (inod>=0) fails the attachment rule, in_geometry is reset to 0. The
node->groups attachment comes from the shared helper
node_attached_element_groups() (group.cc), which scans the active elements
and reads the group of each element from the ELEMENT_GROUP record (NOT from
the ELEMENT record: its first value is the element type).

Rules: -all requires (attached subset of listed) AND (listed subset of
attached) = attached == listed; default/-any/-only require attached subset of
listed (with at least one element). The default/-any/-only inclusion of the
"no element outside the list" rule reproduces the Professional binary
25-10-2023 measured on merge2.dat (a node of group 0 AND group 1 elements is
not merged when only group 0 is listed).

RELATED BEHAVIOR FIX (corpus force12): the edge-load direction of the area()
integral is now always the PHYSICAL side normal computed from the side node
coordinates (area.cc) — with a geometry_point selector the accumulated
geometry() normals are RADIAL and tilt the load away from the true face
normal. Also: elements of pure -empty/-none groups skip the area() integrals
(elem.cc), so a face shared by a material element and an -empty element is
not loaded twice.
