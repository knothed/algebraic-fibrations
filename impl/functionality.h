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
arr2d_fixed reduce_colorings(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos, int num_threads);

// graph.c
arr2d_fixed get_isometries(arr2d_fixed adj);

// legal.c
typedef struct {
    // all colorings for which there are legal orbits
    arr2d_fixed colorings;
    // each coloring can have multiple legal orbits, so this array ontains one row per row in result_cols
    arr2d_var states;
} legal_orbits_result;

// Arguments that the worker threads receive.
typedef struct {
    int n;
    arr2d_fixed legal_states;
    bool* legal_dict;
    arr2d_fixed colorings;
    int start_idx;
    int end_idx;
    legal_orbits_result result;
    bool* stop;
    bool stop_after_first;
    int* num_done;
} orbit_thread_args;

// Return type from find_legal_orbits. When threaded, this lets the user track the calculation progress.
typedef struct {
    // calculation
    int n;
    int num_threads;
    pthread_t* pids;
    orbit_thread_args* args;
    bool* stop;

    // updated via calc_update()
    double progress;
    uint64_t estimated_ms;
    bool finished;

    int num_colorings;
    uint64_t begin_ms;
} legal_orbits_calculation;

arr2d_fixed all_legal_states(arr2d_fixed adj, arr2d_fixed isos);
legal_orbits_calculation find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states, int num_threads, bool force_threaded, bool stop_after_first);

legal_orbits_calculation calc_update(legal_orbits_calculation calc);
legal_orbits_result calc_finish(legal_orbits_calculation calc);

#endif