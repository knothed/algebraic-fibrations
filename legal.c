///
/// Functions related to state and orbit legality.
///

#include <stdio.h>
#include <stdbool.h>
#include "utils.h"

/******* LEGAL STATES *******/

#define MAX_VERTS 32
bool subgraph_connected(arr2d_fixed adj, int sub_size, int vertices[]);
bool is_state_legal(arr2d_fixed adj, int state);

// All legal states between 0 and 2^(n-1).
// To avoid redundancy, we don't return any legal states where vertex (n-1) is set.
// The result is a 2D array where each row is a single integer representing the state as a bitmask.
arr2d_fixed all_legal_states(arr2d_fixed adj) {
    int n = adj.len;
    int max = 1 << (n-1);
    arr2d_fixed result = arr2d_fixed_create_empty(1,max);
    int count = 0;

    for (int i=1; i<max; i++) {
        if (is_state_legal(adj,i)) {
            result.data[count] = i;
            count++;
        }
    }

    result.len = count;
    return result;
}

// Check whether the ascending link and descending link given by the state are both connected and nonempty.
// In other words, check whether the state is legal.
// The graph is given by its adjacency matrix in contiguous form, n successive entries making a row.
// The state is given as a bitmask, with bits 0 to n-1 representing the vertices.
bool is_state_legal(arr2d_fixed adj, int state) {
    int n = adj.len;

    // Read graph columns from state bitmask
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
    return subgraph_connected(adj, asc_size, asc) && subgraph_connected(adj, n-asc_size, desc);
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

bool subgraph_connected(arr2d_fixed adj, int sub_size, int vertices[]) {
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
            if (!visited[i] & get_arrf(adj,vertices[v],vertices[i])) {
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

static inline arr2d_fixed find_legal_orbits_single(int n, int* coloring, arr2d_fixed legal_states, bool* legal, arr2d_fixed append_result);

arr2d_fixed find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states) {
    arr2d_fixed result = arr2d_fixed_create_empty(1, 10);

    // Convert legal_states into "dictionary" for fast checking whether state is contained
    int max_states = 1 << (n-1);
    bool legal_dict[max_states];
    memset(legal_dict,0,max_states*sizeof(bool));
    for (int i=0; i<legal_states.len; i++)
        legal_dict[get_arrf1d(legal_states,i)] = 1;

    // Do each coloring
    for (int i=0; i<colorings.len; i++) {
        bool legal_copy[max_states];
        memcpy(legal_copy,legal_dict,sizeof(legal_copy));
        result = find_legal_orbits_single(n, colorings.data+n*i,legal_states,legal_copy,result);
    }

    return result;
}

// Find all legal orbits from the given legal states in the given coloring.
// The returned list contains a single legal state per legal orbit.
static inline arr2d_fixed find_legal_orbits_single(int n, int* coloring, arr2d_fixed legal_states, bool* legal, arr2d_fixed append_result) {
    arr2d_fixed result = append_result;

    // Convert coloring list into vertex bitmasks
    int color_masks[n];
    int num_cols = 0;
    memset(color_masks,0,n*sizeof(int));
    for (int i=0; i<n; i++) {
        color_masks[coloring[i]] += (1<<i);
        num_cols = MAX(num_cols, coloring[i]+1);
    }

    // Go through all legal states until none are left
    int max_states = 1 << (n-1);
    int idx = 0;
    int remaining = legal_states.len;
    int orbit_size = 1 << num_cols;
    int half_orbit_size = orbit_size >> 1;

    while (remaining >= half_orbit_size) {
        int state = get_arrf1d(legal_states,idx);
        if (!legal[state]) {
            idx++;
            continue;
        }

        // Check if orbit is legal and delete it from the dictionary simultaneously
        bool orbit_legal = true;

        int acted = state;
        int binary = 0;
        for (int c=0; c<orbit_size; c++) {
            if (acted < max_states) { // I/O-symmetry of state; the negative state is in the same orbit
                // Remove from dictionary
                if (legal[acted]) {
                    remaining--;
                    legal[acted] = false;
                } else {
                    orbit_legal = false;
                }
            }

            // Act via gray code
            if (!(c & 1)) {
                binary ^= 1;
                acted ^= color_masks[0];
            } else {
                int y = binary & (~(binary-1));
                binary ^= (y << 1);
                acted ^= color_masks[log2_int(y)+1];
            }
        }

        if (orbit_legal)
            result = append_arrf_single(result, state);
    }

    return result;
}


/*
// Alternative algorithm: reduce all states into a single coset, then sort and analyze this list.
// We need a very quick (radix) sort for this to be efficient!
static inline arr2d_fixed find_legal_orbits_alternative(int n, int* coloring, arr2d_fixed legal_states, arr2d_fixed append_result) {
    arr2d_fixed result = append_result;

    // Convert coloring list into vertex bitmasks
    int color_masks[n];
    int highest_bits[n]; // position of the highest bit of each color
    memset(color_masks,0,n*sizeof(int));
    memset(highest_bits,~0,n*sizeof(int)); // init to -1
    int num_cols = 0;
    for (int i=0; i<n; i++) {
        color_masks[coloring[i]] += (1<<i);
        highest_bits[coloring[i]] = MAX(highest_bits[coloring[i]], i);
        num_cols = MAX(num_cols, coloring[i]+1);
    }

    // Reduce all legal states into one equivalence class by setting all the highest bits to zero
    int states[legal_states.len];

    for (int i=0; i<legal_states.len; i++) {
        int state = legal_states.data[i];
        int elem = 0;
        for (int c=0; c<num_cols; c++) {
            if ((state >> highest_bits[c]) & 1)
                elem |= color_masks[c];
        }
        states[i] = state ^ elem;
    }

    int_radix_sort(states,legal_states.len); // qsort is way too slow
    int half_orbit_size = 1 << (num_cols-1);

    int curr_state = -1;
    int curr_state_begin = 0;
    for (int i=0; i<legal_states.len; i++) {
        if (states[i] != curr_state) {
            if (i-curr_state_begin == half_orbit_size)
                result = append_arrf_single(result, curr_state);
            curr_state = states[i];
            curr_state_begin = i;
        }
    }
    if (legal_states.len-curr_state_begin == half_orbit_size)
        result = append_arrf_single(result, curr_state);

    return result;
}
*/