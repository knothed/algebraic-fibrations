///
/// Functions related to creation of colorings.
///

#include <stdio.h>
#include <stdbool.h>
#include "functionality.h"
#include "sort_r.h"
#include "utils.h"

/******** UPPER BOUND ON NUMBER OF COLORS ********/

// Deduce an upper bound on the number of colors for a coloring from cliques in the graph and the shape of the legal states.
// legal_states only contains non-redundant legal states from 0 to 2^(n-1).
// cliques only contains cliques of size 2 or more.
int num_colors_upper_bound(int n, arr2d_var cliques, arr2d_fixed legal_states) {
    // 1-clique check (number of legal states)
    int upper_bound = log2_int(legal_states.len) + 1;
    int p = upper_bound;

    for (int i=0; i<cliques.len; i++) {
        int size = size_arrv(cliques,i);
        int max = 1<<(size-1);
        int counts[max];
        memset(counts,0,max*sizeof(int));

        // Count how often (all combinations of the bits at the positions of the vertices in the clique) occur in the legal states
        for (int j=0; j<legal_states.len; j++) {
            int bits = 0;
            for (int b=0; b<size; b++) {
                bits += ((get_arrf1d(legal_states,j) >> get_arrv(cliques,i,b)) & 1) * (1<<b); // check if bit b is set in legal_states[j]
            }

            if (bits >= max)
                bits = 2*max - bits - 1;
            counts[bits]++;
        }

        // Evaluate
        int min = legal_states.len;
        for (int c=0; c<max; c++)
            min = MIN(min, counts[c]);

        upper_bound = MIN(upper_bound, log2_int(min << size));
    }

    // if (upper_bound < p)
    //    printf("Reduced from %d to %d!\n", p, upper_bound);

    return upper_bound;
}

/******** FIND ALL COLORINGS ********/

// Greedily convert the list of cliques into a vertex partition such that every partition set is a clique.
// The given cliques must be sorted by length in a descending order.
// Attention: for performance reasons when coloring, just return the first clique, the other vertices are returned as singletons.
arr2d_var cliquewise_vertex_partition(int n, arr2d_var cliques) {
    arr2d_var partition = arr2d_var_create_empty(n,n);
    int current_count = 0;

    for (int i=0; i<1; i++) { // just get a single largest clique
    //for (int i=0; i<cliques.len; i++) {
        int clique_size = size_arrv(cliques,i);
        int clique_overlaps = false;

        if (clique_size > n-current_count) {
            goto out; // assumes the cliques are sorted by length descendingly
        }

        // check if cliques[i] overlaps with partitions[j]
        for (int j=0; j<partition.len; j++) {
            int partition_size = size_arrv(partition,j);
            // probably inefficient, but this is whole procedure is once done once per graph
            for (int k=0; k<clique_size; k++) {
                for (int l=0; l<partition_size; l++) {
                    if (get_arrv(cliques,i,k) == get_arrv(partition,j,l)) {
                        clique_overlaps = true;
                        goto end;
                    }
                }
            }
        }

        end:
        if (!clique_overlaps) {
            int start_idx = (i==0) ? 0 : cliques.end_indices[i-1];
            partition = append_arrv(partition,cliques.data+start_idx,clique_size);
            current_count += clique_size;
        }
    }
    out: if (0) {}

    // Add remaining single vertices
    int v = 0;
    while (current_count < n) {
        start:
        for (int i=0; i<current_count; i++) {
            if (partition.data[i] == v) {
                v++;
                goto start;
            }
        }
        partition = append_arrv_single(partition,v);
        current_count++;
        v++;
    }

    return partition;
}

arr2d_fixed find_all_colorings_impl(arr2d_fixed adj, int num_cols, int used_cols, arr2d_var partition, arr2d_fixed result, int current_coloring[], int level);

// Find all graph colorings with num_cols colors, using the given cliquewise vertex partition for more efficiency.
// Because the colors of one largest clique are fixed, the returned colorings are all pairwise non-equivalent under color relabeling.
// The given cliques must be sorted by length in a descending order.
// Precondition: num_cols <= 32.
arr2d_fixed find_all_colorings(arr2d_fixed adj, int num_cols, arr2d_var partition) {
    int n = adj.len;
    arr2d_fixed result = arr2d_fixed_create_empty(n, 100);
    int current_col[n];
    for (int i=0; i<n; i++) current_col[i] = -1;
    return find_all_colorings_impl(adj, num_cols, 0, partition, result, current_col, 0);
}

arr2d_fixed find_all_colorings_impl(arr2d_fixed adj, int num_cols, int used_cols, arr2d_var partition, arr2d_fixed result, int current_coloring[], int level) {
    int n = adj.len;

    // Add coloring
    if (level == partition.len)
        return append_arrf(result, current_coloring);

    int clique_size = size_arrv(partition,level);
    int remaining = n-partition.end_indices[level];

    // very first clique
    if (level == 0) {
        if (clique_size > num_cols || num_cols > n)
            return result;

        for (int i=0; i<size_arrv(partition,0); i++)
            current_coloring[get_arrv(partition,0,i)] = i;
        return find_all_colorings_impl(adj, num_cols, clique_size, partition, result, current_coloring, level+1);
    }

    // now: level > 0. Invariants:
    // clique_size <= used_cols <= num_cols <= n
    int min_new_cols = MAX(0, num_cols-used_cols-remaining);
    int max_new_cols = MIN(clique_size, num_cols-used_cols);

    for (int new_cols=min_new_cols; new_cols<=max_new_cols; new_cols++) {
        // choose new_cols vertices which get a new color
        int _new_col_verts_local[new_cols * choose(clique_size, new_cols)];
        arr2d_fixed new_col_verts = do_choose(clique_size, new_cols, _new_col_verts_local);

        // go through all options
        for (int i=0; i<new_col_verts.len; i++) {
            int used_indices[new_cols];

            // apply subset to coloring
            for (int j=0; j<new_cols; j++) {
                int idx = get_arrf(new_col_verts,i,j);
                used_indices[j]=idx;
                int v = get_arrv(partition,level,idx);
                current_coloring[v] = used_cols+j;
            }

            // find the remaining vertices
            // we use that `choose` always returns the k numbers in a sorted fashion
            int rem_count = clique_size-new_cols;
            int remaining_indices[rem_count];
            int c1=0;
            int c2=0;
            for (int i=0; i<clique_size; i++) {
                if (c1 < new_cols && used_indices[c1] == i) {
                    c1++;
                } else {
                    remaining_indices[c2] = i;
                    c2++;
                }
            }

            // store the forbidden colors for each remaining vertex in a bitmask (this means that num colors <= 32)
            int forbidden[rem_count];
            memset(forbidden,0,rem_count*sizeof(int));
            for (int r=0; r<rem_count; r++) {
                int v = get_arrv(partition,level,remaining_indices[r]);
                for (int k=0; k<n; k++) {
                    if (get_arrf(adj,v,k) & (current_coloring[k] >= 0))
                        forbidden[r] |= (1<<current_coloring[k]);
                }
            }

            // choose and distribute clique_size-new_cols used colors to the remaining vertices
            int _remaining_vert_cols_local[rem_count * ordered_choose(used_cols, rem_count)];
            arr2d_fixed remaining_vert_cols = do_ordered_choose(used_cols, rem_count, _remaining_vert_cols_local);

            // go through all options
            for (int i=0; i<remaining_vert_cols.len; i++) {
                // legality check
                bool valid = true;
                for (int r=0; r<rem_count; r++) {
                    int c = get_arrf(remaining_vert_cols,i,r);
                    if ((forbidden[r] >> c) & 1) {
                        valid = false;
                        goto out;
                    }
                }
                out:
                if (!valid)
                    continue;

                int new_col[n];
                memcpy(new_col,current_coloring,n*sizeof(int));

                // apply choosing to coloring
                for (int j=0; j<rem_count; j++) {
                    int v = get_arrv(partition,level,remaining_indices[j]);
                    new_col[v] = get_arrf(remaining_vert_cols,i,j);
                }

                result = find_all_colorings_impl(adj, num_cols, used_cols+new_cols, partition, result, new_col, level+1);
            }

            // undo coloring
            for (int j=0; j<new_cols; j++) {
                int idx = get_arrf(new_col_verts,i,j);
                int v = get_arrv(partition,level,idx);
                current_coloring[v] = -1;
            }
        }
    }

    return result;
}


/******** REDUCE COLORINGS BY ISOMETRIES ********/

static inline void make_canonical_form(int n, int* coloring, int num_cols, arr2d_fixed isos);
int cmp_lexicographic(const void* a, const void* b, void* arg);

// Reduce the array of colorings up to color swapping and graph isomorphism.
arr2d_fixed reduce_colorings(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos) {
    // If there is only the identity iso, nothing can be reduced because of the way we generated the colorings
    if (isos.len <= 1) {
        arr2d_fixed copy = arr2d_fixed_create_empty(cols.row_len, cols.len);
        copy = append_arrf_multiple(copy, cols);
        return copy;
    }

    // 1. Bring all colorings into a canonical form: the form obtained by color swapping and graph isos which is lexicographically lowest.
    int64_t a = millis();
    for (int i=0; i<cols.len; i++)
        make_canonical_form(n, cols.data+n*i, num_colors, isos);

    int64_t b = millis();
    sort_r(cols.data, cols.len, n*sizeof(int), cmp_lexicographic, &n);
    printf("can. form took %lld ms, sorting took %lld ms\n", b-a, millis()-b);

    // 2. Remove multiples from the colorings list
    arr2d_fixed result = arr2d_fixed_create_empty(n, cols.len/isos.len);
    int last_new_idx = -1;
    for (int i=0; i<cols.len; i++) {
        if (last_new_idx < 0 || cmp_lexicographic(cols.data+n*i, cols.data+n*last_new_idx, &n)) {
            result = append_arrf(result, cols.data+n*i);
            last_new_idx = i;
        }
    }

    return result;
}

int cmp_lexicographic(const void* a, const void* b, void* arg) {
    int n = *(int*)arg;
    for (int j=0; j<n; j++) {
        int cmp = *(((int*)a)+j) - *(((int*)b)+j);
        if (cmp) return cmp;
    }
    return 0;
}

// Bring the coloring into a canonical form: the form obtained by color swapping and graph isos which is lexicographically lowest.
// This is done in-place.
static inline void make_canonical_form(int n, int* coloring, int num_cols, arr2d_fixed isos) {
    int best[n]; // current lexicographically lowest coloring

    // For each iso, perform color swaps to get the lexicographically lowest coloring for this iso and compare with the current best
    for (int i=0; i<isos.len; i++) {
        int current_col = 0;
        int dict[num_cols]; // mapping old colors -> new colors
        memset(dict,~0,num_cols*sizeof(int)); // init to -1

        int became_better_at = 100;
        bool is_better = (i==0);

        // bring into lexicographically lowest form
        for (int j=0; j<n; j++) {
            int v = get_arrf(isos,i,j);
            int old_col = coloring[v];
            int new_col = dict[old_col];
            if (new_col < 0) {
                dict[old_col] = current_col; // do color swapping assignment
                new_col = current_col;
                current_col++;
            }

            // check wheter current coloring is worse, better (or same) than the current best
            if (!is_better && new_col > best[j]) break;
            if (!is_better && new_col < best[j]) is_better = true;

            // update if better
            if (is_better) best[j] = new_col;
        }
    }

    memcpy(coloring,best,n*sizeof(int));
}