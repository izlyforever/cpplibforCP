<head>
	<script type="text/x-mathjax-config">
		MathJax.Hub.Config({
		  tex2jax: {
			skipTags: ['script', 'noscript', 'style', 'textarea', 'pre'],
			inlineMath: [['$','$']],
			processEscapes: true
		  }
		});
	</script>
	<script type="text/javascript" async
	  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/latest.js?config=TeX-MML-AM_CHTML">
	</script>
</head>

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
- spfa
- Scc: Strongly Connected Components
- twoSAT
- cutVertex
- CutEdge


## Flow

- Dinic: S-T max-Flow $O(n^2 m)$
- HLPP: S-T max-Flow $O(n^2 \sqrt{m})$
- Stoer-Wagner: Global minimum cut of undirected graph($O(n^3)$ implement)
- Flow: minimum cost maximum flow