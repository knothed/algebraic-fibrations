#include <stdio.h>
#include <stdbool.h>

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