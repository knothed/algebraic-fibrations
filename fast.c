#include <stdio.h>
#include <stdbool.h>

/******* ISOMORPHIC COLORINGS *******/

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

/******* GRAPH CONNECTEDNESS *******/

#define MAX_VERTS 20
bool subgraph_connected(int n, int* adj_matrix, int sub_size, int vertices[]);

// Check whether the ascending link and descending link given by the state are both connected and nonempty.
// In other words, check whether the state is legal.
// The graph is given by its adjacency matrix in contiguous form, n successive entries making a row.
// The state is given as a bitmap, with bits 0 to n-1 representing the vertices.
bool is_state_legal(int n, int* adj_matrix, int state) {
    // Read graph columns from state bitset
    int asc_size = 0;
    int asc[n];
    int desc[n];
    for (int k=0; k<n; k++) {
        if ((state >> k) & 1) {
            asc[asc_size] = k;
            asc_size += 1;
        } else {
            desc[k-asc_size] = k;
        }
    }

    if (asc_size == 0 || asc_size == n)
        return false; // subgraph is empty

    // Check connectedness
    return subgraph_connected(n, adj_matrix, asc_size, asc) && subgraph_connected(n, adj_matrix, n-asc_size, desc);
}

typedef struct {
    int queue[MAX_VERTS];
    int front; // = -1
    int rear; // = -1
} bfs_queue;

bool queue_empty(bfs_queue queue) {
    return queue.front == -1 || queue.front > queue.rear;
}

bfs_queue queue_insert(bfs_queue queue, int v) {
    bfs_queue new = queue;

    if (new.front == -1)
        new.front = 0;
    new.rear++;
    new.queue[new.rear] = v;
    return new;
}

bfs_queue queue_delete(bfs_queue queue, int* v) {
    bfs_queue new = queue;

    *v = new.queue[new.front];
    new.front++;
    return new;
}

bool subgraph_connected(int n, int* adj, int sub_size, int vertices[]) {
    bool visited[sub_size];
    memset(visited,0,sub_size*sizeof(bool));

    int v = 0;
    bfs_queue queue = {.front = -1, .rear = -1};
    queue = queue_insert(queue, v); // we label the vertices (0, ..., sub_size-1) and then translate via `vertices`

    // Do BFS
    while (!queue_empty(queue)) {
        queue = queue_delete(queue, &v);
        visited[v] = 1;
        for (int i=0; i<sub_size; i++) { // add adjacent unvisited vertices to queue
            if (!visited[i] && adj[n*vertices[v]+vertices[i]]) {
                visited[i] = 1;
                queue = queue_insert(queue, i);
            }
        }
    }

    // Check if all vertices have been visited
    for (int i=1; i<sub_size; i++) {
        if (!visited[i])
            return false;
    }

    return true;
}


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