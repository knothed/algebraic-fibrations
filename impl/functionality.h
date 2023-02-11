///
/// Exposes functionality of coloring.c, graph.c and legal.c for use by fibering.c
///

#ifndef FUNCTIONALITY_H
#define FUNCTIONALITY_H

#include <stdio.h>
#include <stdbool.h>
#include "utils.h"

// coloring.c
int num_colors_upper_bound(int n, arr2d_var cliques, arr2d_fixed legal_states);
arr2d_var cliquewise_vertex_partition(int n, arr2d_var cliques);
arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partition);
arr2d_fixed reduce_colorings(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos);

// graph.c
arr2d_fixed get_isometries(arr2d_fixed adj);

// legal.c
typedef struct {
    // all colorings for which there are legal orbits
    arr2d_fixed colorings;
    // each coloring can have multiple legal orbits, so this array ontains one row per row in result_cols
    arr2d_var states;
} legal_orbits_result;

arr2d_fixed all_legal_states(arr2d_fixed adj, arr2d_fixed isos);
legal_orbits_result find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states, int num_threads, bool stop_after_first);

#endif