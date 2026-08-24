# NTHU VLSI Physical Design Automation HWs

Implementations for the 114 fall semester NTHU VLSI Physical Design Automation
homework assignments.

The assignments follow a simplified physical-design flow:

1. Partition a circuit to reduce the number of cut nets.
2. Legalize and optimize the placement of standard cells.
3. Route nets on a grid while controlling congestion.

## HW2: Min-cut Partitioning

HW2 partitions the cells of a netlist into either two or four groups. The
objective is to minimize the number of cut nets while keeping the total cell
size of each group within the required balance range.

> **Key keywords:** Multilevel coarsening, FM partitioning, gain bucket, uncoarsening

The implementation uses a multilevel flow:

- Coarsen the netlist by merging cells connected by nets.
- Partition the coarse graph with an FM-style gain-bucket refinement.
- Uncoarsen the graph and refine the partition at each finer level.
- Build a four-way partition through repeated two-way partitioning.


## HW3: Detailed Placement

HW3 performs detailed placement on a valid global placement. It reads the
technology and library information from LEF and the design from DEF, then
moves movable standard cells onto legal rows and sites without overlap. The
main optimization objective is reduced total wirelength.

> **Key keywords:**  global swap, vertical swap, local reordering

The placement flow includes:

- LEF/DEF parsing and conversion to database units.
- Row assignment while preserving fixed components.
- Global swaps and vertical swaps to search nearby rows.
- Local reordering of small groups of cells.
- Single-segment optimization of cell positions within a row.
- Wirelength evaluation using net bounding boxes.


## HW4: Global Routing

HW4 routes two-pin nets on a rectangular routing grid. It tries to minimize
overflow first and wirelength second, subject to the horizontal and vertical
edge capacities supplied in the input.

> **Key keywords:** congestion analysis, rip-up and reroute, A* search, overflow minimization


The routing flow includes:

- Initial L-shaped routes for each two-pin net.
- Congestion and overflow evaluation from edge usage and capacity.
- Prefix sums for fast horizontal and vertical route-cost queries.
- Rip-up and reroute for nets using congested edges.
- A* search for alternative routes.
- Parameter exploration over multiple routing-cost settings.
