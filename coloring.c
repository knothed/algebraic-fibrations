///
/// Functions related to creation of colorings.
///

#include <stdio.h>
#include <stdbool.h>
#include "utils.h"

/******** UPPER BOUND ON NUMBER OF COLORS ********/

int log2_int(int a);

// Deduce an upper bound on the number of colors for a coloring from cliques in the graph and the shape of the legal states.
// legal_states only contains non-redundant legal states from 0 to 2^(n-1).
// cliques only contains cliques of size 2 or more.
// cliques_start_indices must contain cliques_count+1 entries, the first being 0 and the last being the total size of the continuous cliques array.
int max_possible_colors(int n, arr_of_arrs cliques, int* legal_states, int legal_count) {
    // 1-clique check (number of legal states)
    int upper_bound = log2_int(legal_count) + 1;
    int p = upper_bound;

    for (int i=0; i<cliques.len; i++) {
        int size = arr_size(cliques,i);
        int max = 1<<(size-1);
        int counts[max];
        memset(counts,0,max*sizeof(int));

        // Count how often (all combinations of the bits at the positions of the vertices in the clique) occur in the legal states
        for (int j=0; j<legal_count; j++) {
            int bits = 0;
            for (int b=0; b<size; b++) {
                bits += ((legal_states[j] >> arr_get(cliques,i,b)) & 1) * (1<<b); // check if bit b is set in legal_states[j]
            }

            if (bits >= max)
                bits = 2*max - bits - 1;
            counts[bits]++;
        }

        // Evaluate
        int min = legal_count;
        for (int c=0; c<max; c++)
            min = MIN(min, counts[c]);

        upper_bound = MIN(upper_bound, log2_int(min << size));
    }

    if (upper_bound < p)
        printf("Reduced from %d to %d!\n", p, upper_bound);

    return upper_bound;
}

int log2_int(int a) {
    if (a <= 0) return -1;
    int r = 0;
    int b = a;
    while (b >>= 1) r++;
    return r;
}

/******** FIND ALL COLORINGS ********/

// Greedily convert the list of cliques into a vertex partition such that every partition set is a clique.
// The given cliques must be sorted by length in a descending order.
// Attention: for performance reasons when coloring, just return the first clique, the other vertices are returned as singletons.
arr_of_arrs cliquewise_vertex_partition(int n, arr_of_arrs cliques) {
    arr_of_arrs partition = create_arr_of_arrs(n,n);
    int current_count = 0;

    for (int i=0; i<1; i++) { // just get a single largest clique
    //for (int i=0; i<cliques.len; i++) {
        int clique_size = arr_size(cliques,i);
        int clique_overlaps = false;

        if (clique_size > n-current_count) {
            goto out; // assumes the cliques are sorted by length descendingly
        }

        // check if cliques[i] overlaps with partitions[j]
        for (int j=0; j<partition.len; j++) {
            int partition_size = arr_size(partition,j);
            // probably inefficient, but this is whole procedure is once done once per graph
            for (int k=0; k<clique_size; k++) {
                for (int l=0; l<partition_size; l++) {
                    if (arr_get(cliques,i,k) == arr_get(partition,j,l)) {
                        clique_overlaps = true;
                        goto end;
                    }
                }
            }
        }

        end:
        if (!clique_overlaps) {
            partition = arr_append(partition,cliques.data+cliques.start_indices[i],clique_size);
            current_count += clique_size;
        }
    }
    out:
    if (0) {}

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
        int p[1] = {v};
        partition = arr_append(partition,p,1);
        current_count++;
        v++;
    }

    return partition;
}

int* find_all_colorings_impl(int n, int* adj, int num_cols, int used_cols, arr_of_arrs cliques, int* result, int* result_count, int current_coloring[], int level);

// Find all graph colorings with num_cols colors, using the given cliquewise vertex partition for more efficiency.
// Because the colors of one largest clique are fixed, the returned colorings are all pairwise non-equivalent under color relabeling.
// The given cliques must be sorted by length in a descending order.
// Precondition: num_cols <= 32.
int* find_all_colorings(int n, int* adj, int num_cols, arr_of_arrs partition, int* result_count) {
    int size_guess = 100;
    int* result = (int*)malloc(size_guess*n*sizeof(int));
    *result_count = 0;

    int current_col[n];
    for (int i=0; i<n; i++) current_col[i] = -1;
    return find_all_colorings_impl(n, adj, num_cols, 0, partition, result, result_count, current_col, 0);
}

int* find_all_colorings_impl(int n, int* adj, int num_cols, int used_cols, arr_of_arrs partition, int* result, int* result_count, int current_coloring[], int level) {
    int size_guess = 100; // todo: double capacity every time

    // Add coloring
    if (level == partition.len) {
        memcpy(result+n*(*result_count),current_coloring,n*sizeof(int));
        (*result_count)++;
        if ((*result_count)%size_guess == 0)
            result = realloc(result,((*result_count)+size_guess)*n*sizeof(int));
        return result;
    }

    int clique_size = arr_size(partition,level);
    int remaining = n-partition.start_indices[level+1];

    // very first clique
    if (level == 0) {
        if (clique_size > num_cols || num_cols > n)
            return result;

        for (int i=0; i<arr_size(partition,0); i++)
            current_coloring[arr_get(partition,0,i)] = i;
        return find_all_colorings_impl(n, adj, num_cols, clique_size, partition, result, result_count, current_coloring, level+1);
    }

    // now: level > 0. Invariants:
    // clique_size <= used_cols <= num_cols <= n
    int min_new_cols = MAX(0, num_cols-used_cols-remaining);
    int max_new_cols = MIN(clique_size, num_cols-used_cols);

    for (int new_cols=min_new_cols; new_cols<=max_new_cols; new_cols++) {
        // choose new_cols vertices which get a new color
        int count = choose(clique_size, new_cols);
        int new_col_verts[count*new_cols];
        do_choose(clique_size, new_cols, new_col_verts);

        // go through all options
        for (int i=0; i<count; i++) {
            int used_indices[new_cols];

            // apply subset to coloring
            for (int j=0; j<new_cols; j++) {
                int idx = new_col_verts[i*new_cols+j];
                used_indices[j]=idx;
                int v = arr_get(partition,level,idx);
                current_coloring[v] = used_cols+j;
            }

            // find the remaining vertices
            // we use that `choose` always returns the k numbers in a sorted fashion
            int rem_count = clique_size-new_cols;
            int remaining_indices[rem_count];
            int c1=0;
            int c2=0;
            for (int i=0; i<clique_size; i++) {
                if (used_indices[c1] == i) {
                    c1++;
                } else {
                    remaining_indices[c2] = i;
                    c2++;
                }
            }

            // store the forbidden colors for each remaining vertex in a bitset (this means: max colors = 32)
            int forbidden[rem_count];
            memset(forbidden,0,rem_count*sizeof(int));
            for (int r=0; r<rem_count; r++) {
                int v = arr_get(partition,level,remaining_indices[r]);
                for (int k=0; k<n; k++) {
                    if (adj[n*v+k] && current_coloring[k] >= 0)
                        forbidden[r] |= (1<<current_coloring[k]);
                }
            }

            // choose and distribute clique_size-new_cols used colors to the remaining vertices
            int count = ordered_choose(used_cols, rem_count);
            int remaining_vert_cols[count*rem_count];
            do_ordered_choose(used_cols, rem_count, remaining_vert_cols);

            // go through all options
            for (int i=0; i<count; i++) {
                // legality check
                bool valid = true;
                for (int r=0; r<rem_count; r++) {
                    int c = remaining_vert_cols[i*rem_count+r];
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
                    int v = arr_get(partition,level,remaining_indices[j]);
                    new_col[v] = remaining_vert_cols[i*rem_count+j];
                }

                result = find_all_colorings_impl(n, adj, num_cols, used_cols+new_cols, partition, result, result_count, new_col, level+1);
            }

            // undo coloring
            for (int j=0; j<new_cols; j++) {
                int idx = new_col_verts[i*new_cols+j];
                int v = arr_get(partition,level,idx);
                current_coloring[v] = -1;
            }
        }
    }

    return result;
}


/******** REDUCE COLORINGS BY ISOMETRIES ********/

bool is_color_permutation_iso(int n, int c, int* col1, int* col2, int* f);

// Reduce the array of colorings up to color swapping and graph isomorphism.
// cols and isos are contiguous arrays; each entry (coloring or isomorphism) has n values.
int* kill_permutations_and_isos(int n, int num_colors,
                                int* cols, int cols_count,
                                int* isos, int isos_count,
                                int* result_count) {
    int num_result = 0;
    int size_guess = 5;
    int* result = (int*)malloc(size_guess*n*sizeof(int));

    for (int c=0; c<cols_count; c++) {
        bool is_new = 1;
        for (int c2=0; c2<num_result; c2++) {
            for (int f=0; f<isos_count; f++) {
                if (is_color_permutation_iso(n, num_colors, cols+n*c, result+n*c2, isos+n*f)) {
                    is_new = 0;
                    goto out;
                }
            }
        }
        out:
        if (is_new) {
            memcpy(result+n*num_result, cols+n*c, n*sizeof(int));
            num_result++;
            if (num_result%size_guess == 0)
                result = realloc(result, (num_result+size_guess)*n*sizeof(int));
        }
    }

    *result_count = num_result;
    return result;
}

// check whether c2 is some color permutation of c1 under the automorphism f.
bool is_color_permutation_iso(int n, int num_cols ,int* col1, int* col2, int* f) {
    int swaps[num_cols];
    memset(swaps,0,sizeof(int)*num_cols);
    for (int j=0; j<n; j++) {
        int c2 = col2[j];
        if (swaps[c2]) {
            if (swaps[c2] != col1[f[j]]+1)
                return false;
        } else {
            swaps[c2] = col1[f[j]]+1;
        }
    }
    return true;
}