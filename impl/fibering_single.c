///
/// Graph fibering driver code that combines lower-level functionality to perform a fibering check for a single graph.
///

#include <stdio.h>
#include <stdbool.h>
#include "functionality.h"
#include "utils.h"

legal_orbits_result do_orbit_search(int n, arr2d_fixed legal_states, arr2d_fixed colorings, bool verbose, char* text, int num_threads, bool single_orbit, uint64_t begin_time);

// Find all (or just one) legal orbit(s) for the given graph.
// Consider all colorings with between min_cols and max_cols colors. Split up the parallelizable work into num_thread threads.
// When `total_progress_bar`, first all colorings are calculated and then all orbits are searched. Otherwise, the orbits are searched per number of colors, each giving a new progress bar.
legal_orbits_result graph_fiberings(arr2d_fixed adj, arr2d_var cliques, int min_cols, int max_cols, bool verbose, bool total_progress_bar, int num_threads, bool single_orbit) {
    int n = adj.row_len;

    int64_t begin_time;
    char text[128];
    if (verbose) {
        printf("Preparations ... ");
        fflush(stdout);
        begin_time = millis();
    }

    // Preparations
    //int64_t t0 = millis();
    arr2d_fixed isos = get_isometries(adj);
    //int64_t t1 = millis();
    arr2d_fixed legal_states = all_legal_states(adj,isos);
    //int64_t t2 = millis();
    arr2d_var partitions = cliquewise_vertex_partition(n, cliques);
    //int64_t t3 = millis();

    int cmin = size_arrv(partitions,0); // sage's g.chromatic_number() is too slow, so we start lower
    int cmax = num_colors_upper_bound(n, cliques, legal_states);
    int64_t t4 = millis();
    cmin = MAX(cmin, min_cols);
    if (max_cols > 0) cmax = MIN(cmax, max_cols);

    //printf("\isos: %lld, legal: %lld, parts: %lld, max_cols: %lld\n",t1-t0,t2-t1,t3-t2,t4-t3);

    if (verbose) {
        char* legal_len = pretty_int(legal_states.len);
        char* isos_len = pretty_int(isos.len);
        char* duration = pretty_ms(millis()-begin_time, true);
        printf("%s legal states, %s isos; #colors <= %d (took %s).\n", legal_len, isos_len, cmax, duration);
        free(legal_len);
        free(isos_len);
        free(duration);
    }

    legal_orbits_result all_orbits = { .colorings = arr2d_fixed_create_empty(n, 10), .states = arr2d_var_create_empty(20, 10) };

    // OPTION 1: First collect all colorings, then search for legal orbits

    if (total_progress_bar) {
        arr2d_fixed all_reduced = arr2d_fixed_create_empty(n, 10);

        begin_time = millis();
        for (int c=cmin; c<=cmax; c++) {
            if (verbose) {
                printf("\rSearching %d-colorings...", c);
                fflush(stdout);
            }

            // Find colorings
            arr2d_fixed cols = find_all_colorings(adj, c, partitions);
            arr2d_fixed reduced = reduce_colorings(n, c, cols, isos, num_threads);
            all_reduced = append_arrf_multiple(all_reduced, reduced);
            free_arrf(cols);
            free_arrf(reduced);
        }

        if (verbose) {
            char* num_cols = pretty_int(all_reduced.len);
            snprintf(text,128,"\rTesting %s colorings: ", num_cols);
            free(num_cols);
            printf(text);
            fflush(stdout);
        }

        // Do orbit search
        all_orbits = do_orbit_search(n, legal_states, all_reduced, verbose, text, num_threads, single_orbit, begin_time);
        free_arrf(all_reduced);
    }

    // OPTION 2: Search for legal orbits per #colors consecutively

    else {
        for (int c=cmin; c<=cmax; c++) {
            if (verbose) {
                printf("\rTesting %d colors...", c);
                fflush(stdout);
                begin_time = millis();
            }

            // Find colorings
            arr2d_fixed cols = find_all_colorings(adj, c, partitions);
            arr2d_fixed reduced = reduce_colorings(n, c, cols, isos, num_threads);
            free_arrf(cols);

            if (reduced.len == 0) {
                free_arrf(reduced);
                continue;
            }

            if (verbose) {
                char* reduced_len = pretty_int(reduced.len);
                snprintf(text,128,"\rTesting %s %d-colorings: ", reduced_len, c);
                free(reduced_len);
                printf(text);
                fflush(stdout);
            }

            // Do orbit search
            legal_orbits_result orbits = do_orbit_search(n, legal_states, reduced, verbose, text, num_threads, single_orbit, begin_time);
            bool found_orbit = orbits.colorings.len > 0;
            all_orbits.colorings = append_arrf_multiple(all_orbits.colorings, orbits.colorings);
            all_orbits.states = append_arrv_multiple(all_orbits.states, orbits.states);
            free_arrf(reduced, orbits.colorings);
            free_arrv(orbits.states);

            if (found_orbit && single_orbit)
                break;
        }
    }

    free_arrf(isos, legal_states);
    free_arrv(partitions);

    return all_orbits;
}

// Perform an orbit search with an optional progress indicator.
legal_orbits_result do_orbit_search(int n, arr2d_fixed legal_states, arr2d_fixed colorings, bool verbose, char* text, int num_threads, bool single_orbit, uint64_t begin_time) {
    bool progress_indicator = false;
    if (verbose) {
        uint64_t load_per_thread = (uint64_t)legal_states.len * (uint64_t)colorings.len / (uint64_t)MAX(1,num_threads);
        progress_indicator = load_per_thread > 500000000; // ~ 10s
    }

    legal_orbits_calculation orbit_calc = find_legal_orbits(n, colorings, legal_states, num_threads, progress_indicator, single_orbit);

    // Show progress indicator
    if (progress_indicator) {
        do {
            delay(1000);
            orbit_calc = calc_update(orbit_calc);
            print_progress(text, orbit_calc.progress, orbit_calc.estimated_ms);
        } while (!orbit_calc.finished);
    }

    legal_orbits_result orbits = calc_finish(orbit_calc);
    bool found_orbit = orbits.colorings.len > 0;

    if (verbose) {
        char* num_states = pretty_int(total_len_arrv(orbits.states));
        char* num_colorings = pretty_int(orbits.colorings.len);
        char* duration = pretty_ms(millis()-begin_time,true);
        printf(text);
        if (found_orbit)
            printf("found %s legal orbits on %s colorings (took %s).\n", num_states, num_colorings, duration);
        else
            printf("no orbit found (took %s).\n", duration);
        free(num_states);
        free(num_colorings);
        free(duration);
    }

    return orbits;
}