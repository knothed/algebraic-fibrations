///
/// Graph fibering driver code that combines lower-level functionality in the correct way.
///

#include <stdio.h>
#include <stdbool.h>
#include "functionality.h"
#include "utils.h"

// Find all (or just one) legal orbit(s) for the given graph.
// Consider all colorings with up to max_cols colors. Split orbit search work into num_thread threads.
legal_orbits_result graph_fiberings(arr2d_fixed adj, arr2d_var cliques, int min_cols, int max_cols, bool verbose, int num_threads, bool single_orbit) {
    int n = adj.row_len;

    // Preparations
    if (verbose) {
        printf("Preparations ... ");
        fflush(stdout);
    }
    int64_t now = millis();

    arr2d_fixed isos = get_isometries(adj);
    arr2d_fixed legal_states = all_legal_states(adj,isos);
    arr2d_var partitions = cliquewise_vertex_partition(n, cliques);

    int cmin = size_arrv(partitions,0); // sage's g.chromatic_number() is too slow, so we start lower
    int cmax = num_colors_upper_bound(n, cliques, legal_states);
    cmin = MAX(cmin, min_cols);
    if (max_cols > 0) cmax = MIN(cmax, max_cols);

    if (verbose)
        printf("took %lld ms: %d legal states, %d isos; #colors <= %d\n", millis()-now, legal_states.len, isos.len, cmax);

    // Colorings: all colorings of specific #colors
    legal_orbits_result all_orbits = { .colorings = arr2d_fixed_create_empty(n, 10), .states = arr2d_var_create_empty(20, 10) };

    for (int c=cmin; c<=cmax; c++) {
        if (verbose) {
            printf("Testing %d colors... ", c);
            fflush(stdout);
        }

        // Find colorings
        now = millis();
        arr2d_fixed cols = find_all_colorings(adj, c, partitions);
        arr2d_fixed reduced = reduce_colorings(n, c, cols, isos, num_threads);
        free_arrf(cols);

        if (reduced.len == 0) {
            free_arrf(reduced);
            if (verbose) printf("\n");
            continue;
        }

        if (verbose) {
            printf("found %d unique colorings (in %lld ms)... ", reduced.len, millis()-now);
            fflush(stdout);
        } else if (verbose) {
            printf("\n");
        }

        // Find legal orbits
        now = millis();
        legal_orbits_result orbits = find_legal_orbits(n, reduced, legal_states, num_threads, single_orbit);
        all_orbits.colorings = append_arrf_multiple(all_orbits.colorings, orbits.colorings);
        all_orbits.states = append_arrv_multiple(all_orbits.states, orbits.states);

        // TODO: postprocess orbits: reduce colorings inside orbits if they weren't reduced beforehand!

        bool found_orbit = orbits.colorings.len > 0;
        if (verbose) {
            if (found_orbit)
                printf("found %d legal orbits on %d colorings (in %lld ms)!\n", total_len_arrv(orbits.states), orbits.colorings.len, millis()-now);
            else
                printf("found no legal orbit (in %lld ms).\n", millis()-now);
        }

        free_arrf(reduced, orbits.colorings);
        free_arrv(orbits.states);

        if (found_orbit && single_orbit)
            break;
    }

    free_arrf(isos, legal_states);
    free_arrv(partitions);

    return all_orbits;
}