# algebraic-fibrations

This is a tool I wrote for my masters thesis. It works together with Sage to provide functionality for checking whether graphs of your choice emit a coloring-induced legal orbit which result in a virtual fibration of the respective right-angled Coxeter group, as introduced by [Jankiewicz, Norin and Wise in 2017](https://arxiv.org/abs/1711.11505).

The code is inspired by the paper of [Italiano, Martelli and Migliorini in 2020](https://arxiv.org/abs/2010.10200) and explained in more detail in [my thesis](todo).

Given any (reasonably small) graph $\Gamma$, this tool searches through all nonisomorphic colorings of $\Gamma$ and checks whether any of these emits a legal orbit.
For maximum performance, this brute-force search is written in C and is multi-threaded. All the relevant functionality is exposed to Python.

### How to use `algebraic-fibrations`:

 1. Download this tool
 2. Start `sage` in this directory
 3. Call `load('main.sage')`.

Now take any graph in Sage, such as:

```sage
g = polytopes.cuboctahedron().graph()
```
Now find one or all legal orbits:

```sage
one_legal_orbit(g)
all_legal_orbits(g)
```

You can customize the search by specifiying the range of colors of the colorings and the number of threads used for parallelizable operations.

More functionality:

 - `all_reduced_colorings(graph,d)` finds all nonisomorphic d-colorings for `graph`.
 - Use `legal_orbits_for_coloring` to check whether a single coloring emits a legal orbit.
 - Use `graph_k_fibers` whether a graph k-fibers for k > 0.

#### Stream Search

`geng` is a tool that generates a list of graphs with certain properties. With `analyze_geng_stream` you can call geng and analyze all graphs for a legal orbit.

You can also analyze any stream of graphs using `analyze_stream`. The graphs must be separated by newlines and be given in graph6 format.
Use `graph6_from_graph` and `graph_from_graph6` to convert between Sage graphs and the graph6 representations.

### Limitations

Graphs are limited to 32 vertices, but in practice, 20 vertices is the absolute maximum that will work in reasonable time. This is because the number of colorings and the number of subsets of the vertices both grow at least exponentially.