///
/// Functions related to state and orbit legality.
///

#include <stdio.h>
#include <stdbool.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/******* LEGAL STATES *******/

#define MAX_VERTS 40
bool subgraph_connected(int n, int* adj_matrix, int sub_size, int vertices[]);
bool is_state_legal(int n, int* adj_matrix, int state);

// All legal states between 0 and 2^(n-1).
// To avoid redundancy, we don't return any legal states where vertex (n-1) is set.
int* all_legal_states(int n, int* adj_matrix, int* result_count) {
    int max = 1 << (n-1);
    int* result = (int*)malloc(max*sizeof(int));
    int count = 0;

    for (int i=1; i<max; i++) {
        if (is_state_legal(n,adj_matrix,i)) {
            result[count] = i;
            count++;
        }
    }

    *result_count = count;
    return result;
}

// Check whether the ascending link and descending link given by the state are both connected and nonempty.
// In other words, check whether the state is legal.
// The graph is given by its adjacency matrix in contiguous form, n successive entries making a row.
// The state is given as a bitset, with bits 0 to n-1 representing the vertices.
bool is_state_legal(int n, int* adj, int state) {
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
    return subgraph_connected(n, adj, asc_size, asc) && subgraph_connected(n, adj, n-asc_size, desc);
}

/******* GRAPH CONNECTEDNESS *******/

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

/******** LEGAL ORBITS ********/

// Find all legal orbits from the given legal states in the given coloring.
// The returned list contains a single legal state per legal orbit.
int* find_legal_orbits(int n, int* coloring, int* legal_states, int num_states, int* result_count) {
    int size_guess = 5;
    int* result = (int*)malloc(size_guess*sizeof(int));
    *result_count = 0;

    // Convert coloring list into vertex bitmasks
    int color_masks[n];
    int num_cols = 0;
    memset(color_masks,0,n*sizeof(int));
    for (int i=0; i<n; i++) {
        color_masks[coloring[i]] += (1<<i);
        num_cols = MAX(num_cols, coloring[i]+1);
    }

    // Convert legal_states into "dictionary" for fast checking whether state is contained
    int max_states = 1 << (n-1);
    int legal[max_states];
    memset(legal,0,max_states*sizeof(int));
    for (int i=0; i<num_states; i++) {
        legal[legal_states[i]] = 1;
    }

    // Go through all legal states until none are left
    int idx = 0;
    int remaining = num_states;
    int orbit_size = 1 << num_cols;

    while (remaining >= (orbit_size >> 1)) {
        int state = legal_states[idx];
        if (!legal[state]) {
            idx++;
            continue;
        }

        // Check if orbit is legal and delete it from the dictionary simultaneously
        bool orbit_legal = true;

        for (int c=0; c<orbit_size; c++) {
            // Act on state
            int acted = state;
            for (int b=0; b<num_cols; b++) {
                acted ^= ((c >> b) & 1) * color_masks[b]; // act if bit b is set in c
            }

            if (acted >= max_states)
                continue; // I/O-symmetry of state; the negative state is in the same orbit

            // Remove from dictionary
            if (legal[acted]) {
                remaining--;
                legal[acted] = false;
            } else {
                orbit_legal = false;
            }
        }

        if (orbit_legal) {
            result[*result_count] = state;
            (*result_count)++;
            if ((*result_count)%size_guess == 0)
                result = realloc(result, ((*result_count)+size_guess)*sizeof(int));
        }
    }

    return result;
}