///
/// Functions related to graph isometries.
///

#include <stdio.h>
#include <stdbool.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/******** GRAPH ISOMETRIES ********/

int* get_isometries_impl(int n, int* adj, int* isometries, int* isos_count, int current_iso[], int level);

// Calculate all isometries of the graph.
// The graph is given by its adjacency matrix in contiguous form, n successive entries making a row.
int* get_isometries(int n, int* adj, int* result_count) {
    int size_guess = 10;
    int* result = (int*)malloc(size_guess*n*sizeof(int));

    int current_iso[n];
    memset(current_iso,0,n*sizeof(int));

    *result_count = 0;
    result = get_isometries_impl(n, adj, result, result_count, current_iso, 0);
    return result;
}

int* get_isometries_impl(int n, int* adj, int* isometries, int* isos_count, int current_iso[], int level) {
    int size_guess = 10;

    // Add isometry
    if (level == n) {
        memcpy(isometries+n*(*isos_count),current_iso,n*sizeof(int));
        (*isos_count)++;
        if ((*isos_count)%size_guess == 0)
            isometries = realloc(isometries, ((*isos_count)+size_guess)*n*sizeof(int));
        return isometries;
    }

    // insert 'i' somewhere into current_iso
    for (int i=0; i<n; i++) {
        // check whether current_iso already contains i
        bool contains_i = false;
        for (int j=0; j<level; j++) {
            if (current_iso[j] == i) {
                contains_i = true;
                break;
            }
        }

        if (contains_i)
            continue;

        // check whether adding i at position level would preserve all edges
        bool is_ok = true;
        for (int j=0; j<level; j++) {
            if (adj[level*n+j] != adj[i*n+current_iso[j]]) {
                is_ok = false;
                break;
            }
        }

        // add i to current_iso and recurse
        if (is_ok) {
            current_iso[level] = i;
            isometries = get_isometries_impl(n, adj, isometries, isos_count, current_iso, level+1);
        }
    }

    return isometries;
}