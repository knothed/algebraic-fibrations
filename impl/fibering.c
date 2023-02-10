///
/// Graph fibering driver code that combines lower-level functionality in the correct way.
///

#include <stdio.h>
#include <stdbool.h>
#include "functionality.h"
#include "utils.h"

// Find all (or just one) legal orbit(s) for the given graph.
// Consider all colorings with up to max_cols colors. Split orbit search work into num_thread threads.
legal_orbits_result graph_fiberings(arr2d_fixed adj, arr2d_var cliques, int max_cols, bool verbose, int num_threads, bool single_orbit) {
    int n = adj.row_len;

    // Preparations
    if (verbose) {
        printf("Preparations ... ");
        fflush(stdout);
    }
    int64_t now = millis();

    arr2d_fixed legal_states = all_legal_states(adj);
    arr2d_fixed isos = get_isometries(adj);
    arr2d_var partitions = cliquewise_vertex_partition(n, cliques);

    int cmin = 0; // g.chromatic_number() is slow
    int cmax = num_colors_upper_bound(n, cliques, legal_states);
    if (max_cols > 0) cmax = MIN(cmax, max_cols);

    if (verbose)
        printf("took %lld ms: %d legal states, %d isos; #colors <= %d\n", millis()-now, legal_states.len, isos.len, cmax);

    // Colorings: all colorings of specific #colors
    legal_orbits_result all_orbits = { .colorings = arr2d_fixed_create_empty(n, 10), .states = arr2d_var_create_empty(20, 10) };

    // Above REDUCE_THRESHOLD number of colorings, we don't reduce the colorings into isomorpishm classes but directly perform orbit search.
    // Why? The runtime of orbit search is linear in the number of colorings, while the runtime of isomorphism reduction is quadratic.
    // TODO: reduce colorings into cosets -> get O(n log n)
    int REDUCE_THRESHOLD = 20000;

    for (int c=cmin; c<=cmax; c++) {
        if (verbose) {
            printf("Testing %d colors... ", c);
            fflush(stdout);
        }

        // Find colorings
        now = millis();
        arr2d_fixed cols = find_all_colorings(adj, c, partitions);
        if (cols.len < REDUCE_THRESHOLD) {
            arr2d_fixed reduced = kill_permutations_and_isos(n, c, cols, isos);
            free_arrf(cols);
            cols = reduced;
        }

        if (verbose) {
            if (cols.len < REDUCE_THRESHOLD)
                printf("found %d reduced colorings (in %lld ms)... ", cols.len, millis()-now);
            else
                printf("found %d colorings... ", cols.len);
            fflush(stdout);
        }

        // Find legal orbits
        now = millis();
        legal_orbits_result orbits = find_legal_orbits(n, cols, legal_states, num_threads, single_orbit);
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

        free_arrf(cols);
        free_arrf(orbits.colorings);
        free_arrv(orbits.states);

        if (found_orbit && single_orbit)
            break;
    }

    free_arrf(legal_states);
    free_arrf(isos);
    free_arrv(partitions);

    return all_orbits;
}