///
/// Functions related to creation of colorings.
///

#include <stdio.h>
#include <stdbool.h>
#include "utils.h"

/******** UPPER BOUND ON NUMBER OF COLORS ********/

// Deduce an upper bound on the number of colors for a coloring from cliques in the graph and the shape of the legal states.
// legal_states only contains non-redundant legal states from 0 to 2^(n-1).
// cliques only contains cliques of size 2 or more.
// cliques_start_indices must contain cliques_count+1 entries, the first being 0 and the last being the total size of the continuous cliques array.
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

    if (upper_bound < p)
        printf("Reduced from %d to %d!\n", p, upper_bound);

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
            partition = append_arrv(partition,cliques.data+cliques.start_indices[i],clique_size);
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
    int remaining = n-partition.start_indices[level+1];

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

bool precheck(arr2d_var sp1, arr2d_var sp2, int* f);
bool is_color_permutation_iso(int n, int c, int* col1, int* col2, int* f);

typedef struct {
    int idx;
    int shape;
    arr2d_var single_parts; // parts of the partition which occur only once; these are used to eliminate isometries
} col_shape;

// A numeral representation of a partition.
// For n<24, this is injective, i.e. different partitions of n yield different numbers; for n>=24, there are few clashes.
// todo: make this work for n>=24
int shape(int* partition, int parts) {
    int res = 0;
    int pow = 0;

    for (int i=parts-1; i>=0; i--) {
        int v = partition[i];
        if (v>1)
            res += v << pow;
        pow += log2_int(v)+1;
    }

    return res;
}

int cmp_shapes(const void* a, const void* b) {
    return (*(col_shape*)a).shape - (*(col_shape*)b).shape;
}

int cmp_desc(const void* a, const void* b) {
    return *(int*)b - *(int*)a;
}

// Reduce the array of colorings up to color swapping and graph isomorphism.
// cols and isos are contiguous arrays; each entry (coloring or isomorphism) has n values.
arr2d_fixed kill_permutations_and_isos(int n, int num_colors, arr2d_fixed cols, arr2d_fixed isos) {
    if (num_colors >= 24) {
        printf("Error: reducing colorings not yet supported for 24 colors or more");
        exit(1);
    }

    // if there is only the identity iso, nothing can be reduced because of the way we generated the colorings
    if (isos.len <= 1) {
        arr2d_fixed copy = arr2d_fixed_create_empty(cols.row_len, cols.len);
        copy = append_arrf_multiple(copy, cols);
        return copy;
    }

    // 1. categorize colorings by their shape (i.e. partition of n)
    col_shape shapes[cols.len];
    for (int i=0; i<cols.len; i++) {
        // count color occurrences to make partition of num_colors
        int occurrences[num_colors]; // count occurrences of each color
        int vert_occurrences[n]; // count occurrences of the color each vertex has
        memset(occurrences,0,sizeof(int)*num_colors);
        for (int j=0; j<n; j++)
            occurrences[get_arrf(cols,i,j)]++;
        for (int j=0; j<n; j++)
            vert_occurrences[j] = occurrences[get_arrf(cols,i,j)];

        qsort(occurrences,num_colors,sizeof(int),cmp_desc);

        // find parts that occur only once
        arr2d_var single_parts = arr2d_var_create_empty(n,num_colors);
        for (int c=num_colors-1; c>=0; c--) {
            int p = occurrences[c];

            if (c != 0 && p == 1 && p != occurrences[c-1])
                goto ones_exception; // also treat all the colors occurring once as a single part
            if (!(c == 0 || p != occurrences[c-1]) || !(c == num_colors-1 || p != occurrences[c+1]))
                continue;

            ones_exception: if (0) {}
            if (p == 1)
                p = num_colors-c;

            // fill verts with the vertices of the given color
            int verts[p];
            int ix = 0;
            for (int j=0; j<n; j++) {
                if (vert_occurrences[j] == occurrences[c]) {
                    verts[ix] = j;
                    ix++;
                }
            }
            single_parts = append_arrv(single_parts, verts, p);
        }

        shapes[i] = (col_shape) { .idx = i, .shape = shape(occurrences,num_colors), .single_parts = single_parts };
    }

    qsort(shapes,cols.len,sizeof(col_shape),cmp_shapes);

    // 2. only reduce colorings of the same shape
    arr2d_fixed result = arr2d_fixed_create_empty(n, 10);
    arr2d_fixed new_cols = arr2d_fixed_create_empty(n, 5);
    arr2d_fixed new_col_indices = arr2d_fixed_create_empty(1, 5); // index of coloring in shapes
    int prev_shape = -1;

    int prev_shape_idx=0;
    for (int c=0; c<cols.len; c++) {
        // shape finished: copy new colorings of previous shape into result
        if (c > 0 && shapes[c].shape != prev_shape) {
            result = append_arrf_multiple(result,new_cols);
            prev_shape_idx = c;
            new_cols.len = 0;
            new_col_indices.len = 0;
        }
        prev_shape = shapes[c].shape;

        // search for equivalent colorings
        int idx = shapes[c].idx;
        bool is_new = true;

        arr2d_var sp1 = shapes[c].single_parts;
        for (int f=1; f<isos.len; f++) { // skip identity iso
            for (int d=0; d<new_cols.len; d++) {
                arr2d_var sp2 = shapes[get_arrf1d(new_col_indices,d)].single_parts;
                if (!precheck(sp1, sp2, isos.data+n*f))
                    continue;
                if (is_color_permutation_iso(n, num_colors, cols.data+n*idx, new_cols.data+n*d, isos.data+n*f)) {
                    is_new = false;
                    goto out;
                }
            }
        }
        out:
        if (is_new) {
            new_cols = append_arrf(new_cols, cols.data+n*idx);
            new_col_indices = append_arrf_single(new_col_indices,c); // todo: use int* once arrf's use chars
        }
    }

    result = append_arrf_multiple(result,new_cols);
    free_arrf(new_cols);
    return result;
}

// check whether c2 can be a permutation of c1 under f.
bool precheck(arr2d_var sp1, arr2d_var sp2, int* f) {
    for (int s=0; s<sp1.len; s++) {
        int size = size_arrv(sp1,s);
        for (int s1=0; s1<size; s1++) {
            int val = f[get_arrv(sp1,s,s1)];
            bool matched = false;
            for (int s2=0; s2<size; s2++) {
                if (get_arrv(sp2,s,s2) == val) {
                    matched = true;
                    goto matched_out;
                }
            }
            matched_out:
            if (!matched)
                return false;
        }
    }
    return true;
}

// check whether c2 is some color permutation of c1 under the automorphism f.
bool is_color_permutation_iso(int n, int num_cols, int* col1, int* col2, int* f) {
    int swaps[num_cols];
    memset(swaps,0,sizeof(int)*num_cols);
    for (int j=0; j<n; j++) {
        int c1 = col1[j];
        if (swaps[c1]) {
            if (swaps[c1] != col2[f[j]]+1)
                return false;
        } else {
            swaps[c1] = col2[f[j]]+1;
        }
    }
    return true;
}