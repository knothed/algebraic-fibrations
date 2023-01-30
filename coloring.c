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