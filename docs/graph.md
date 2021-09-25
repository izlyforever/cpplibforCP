# graph.hpp

Notation: $n$(numbers of vertex), $m$(numbers of edges)

## Tree

- DfsTour
- EulerTour
- LCA
- Prim: Minimum Spanning Tree
- LiuZhu: Minimum tree diagram
- TopSort
- EulerPath

`link/cut Tree`, `dsu on tree` need case-by-case analysis.


## Shortest Path

- Floyd
- BellmanFord
- Dijkstra
- spfa

## Connectivity

- Scc: Strongly Connected Components
- twoSAT
- cutVertex
- CutEdge


## Flow

- Dinic: S-T max-Flow $O(n^2 m)$
- HLPP: S-T max-Flow $O(n^2 \sqrt{m})$
- Stoer-Wagner: Global minimum cut of undirected graph($O(n^3)$ implement)
- Flow: minimum cost maximum flow


## Mixed

- circle3count: count $(a, b), (b, c), (c, a)$ in undirected graph