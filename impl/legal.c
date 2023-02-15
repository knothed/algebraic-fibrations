///
/// Functions related to state and orbit legality.
///

#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include "functionality.h"
#include "utils.h"

/******* LEGAL STATES *******/

#define MAX_VERTS 32
bool subgraph_connected(arr2d_fixed adj, int sub_size, int vertices[]);
bool is_state_legal(arr2d_fixed adj, int state);

// All legal states between 0 and 2^(n-1).
// To avoid redundancy, we don't return any legal states where vertex (n-1) is set.
// The result is a 2D array where each row is a single integer representing the state as a bitmask.
arr2d_fixed all_legal_states(arr2d_fixed adj, arr2d_fixed isos) {
    int n = adj.len;
    int max = 1 << (n-1);
    arr2d_fixed result = arr2d_fixed_create_empty(1,max/2);

    char* dict = malloc(max*sizeof(char)); // for isometries
    memset(dict,~0,max*sizeof(char));

    for (int state=1; state<max; state++) {
        if (dict[state] >= 0)
            continue;

        bool legal = is_state_legal(adj,state);

        // act on state
        for (int iso=0; iso<isos.len; iso++) {
            int acted = 0;
            for (int i=0; i<n; i++)
                if ((state >> i) & 1)
                    acted += (1 << get_arrf(isos,iso,i));
            if (acted >= max) acted = 2*max-acted-1;
            if (dict[acted] >= 0) continue; // orbit might not be free
            dict[acted] = legal ? 1 : 0;
            if (legal)
                result = append_arrf_single(result, acted);
        }
    }

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

static inline legal_orbits_result find_legal_orbits_single(int n, int* coloring, arr2d_fixed legal_states, bool* legal, legal_orbits_result append_result, bool* stop, bool stop_after_first);

void* orbit_thread_enter(void* arg) {
    orbit_thread_args args = *((orbit_thread_args*)arg);
    legal_orbits_result result = args.result;

    int n = args.n;
    int max_states = 1 << (n-1);
    bool* legal_copy = malloc(max_states*sizeof(bool));
    for (int i=args.start_idx; i<=args.end_idx; i++) {
        if (*args.stop) break;
        memcpy(legal_copy,args.legal_dict,max_states*sizeof(bool));
        result = find_legal_orbits_single(n, args.colorings.data+n*i, args.legal_states, legal_copy, result, args.stop, args.stop_after_first);
        (*args.num_done)++;
    }

    ((orbit_thread_args*)arg)->result = result;
}

// Update `percentage` and `estimated_ms` of the calculation by reviewing the threads' done work.
legal_orbits_calculation calc_update(legal_orbits_calculation calc) {
    int sum = 0;
    for (int i=0; i<calc.num_threads; i++)
        sum += *calc.args[i].num_done;

    calc.finished = sum == calc.num_colorings || *calc.stop;
    calc.progress = (double)sum / (double)calc.num_colorings;

    uint64_t taken = millis()-calc.begin_ms;
    if (calc.progress > 0.001) {
        calc.estimated_ms = (uint64_t) (((double)taken) * (1.0-calc.progress) / calc.progress);
    } else {
        calc.estimated_ms = 0;
    }

    return calc;
}

// Wait for the worker threads to finish and return the result.
// Frees any additional memory that was used by the threads.
legal_orbits_result calc_finish(legal_orbits_calculation calc) {
    for (int i=0; i<calc.num_threads; i++)
        pthread_join(calc.pids[i], 0);

    // Collect result
    int result_len = 0;
    int total_states = 0;
    for (int i=0; i<calc.num_threads; i++) {
        result_len += calc.args[i].result.colorings.len;
        total_states += total_len_arrv(calc.args[i].result.states);
    }
    legal_orbits_result result = { .colorings = arr2d_fixed_create_empty(calc.n, result_len), .states = arr2d_var_create_empty(total_states, result_len) };
    for (int i=0; i<calc.num_threads; i++) {
        result.colorings = append_arrf_multiple(result.colorings, calc.args[i].result.colorings);
        result.states = append_arrv_multiple(result.states, calc.args[i].result.states);
    }

    // Free resources
    for (int i=0; i<calc.num_threads; i++) {
        free_arrf(calc.args[i].result.colorings);
        free_arrv(calc.args[i].result.states);
        free(calc.args[i].num_done);
    }
    free(calc.stop);
    free(calc.pids);
    free(calc.args[0].legal_dict);
    free(calc.args);

    return result;
}

// Find all legal orbits from the given legal states in all of the given colorings.
// If required, split up the work on multiple threads.
// The returned list contains one row per orbit, the first n entries being the coloring and the last entry being one state in the orbit.
legal_orbits_calculation find_legal_orbits(int n, arr2d_fixed colorings, arr2d_fixed legal_states, int num_threads, bool force_threaded, bool stop_after_first) {
    // Convert legal_states into "dictionary" for fast checking whether state is contained
    int max_states = 1 << (n-1);
    bool* legal_dict = malloc(max_states*sizeof(bool));
    memset(legal_dict,0,max_states*sizeof(bool));
    for (int i=0; i<legal_states.len; i++)
        legal_dict[get_arrf1d(legal_states,i)] = 1;

    // Create threads
    num_threads = MAX(1, num_threads);
    pthread_t* pids = malloc(num_threads*sizeof(pthread_t));
    orbit_thread_args* args = malloc(num_threads*sizeof(orbit_thread_args));
    bool* stop = malloc(sizeof(bool));
    *stop = false;

    for (int i=0; i<num_threads; i++) {
        int from = (i*colorings.len)/num_threads;
        int to = ((i+1)*colorings.len)/num_threads - 1;
        int* num_done = malloc(sizeof(int));
        *num_done = 0;
        legal_orbits_result result = { .colorings = arr2d_fixed_create_empty(n, 10), .states = arr2d_var_create_empty(20, 10) };
        args[i] = (orbit_thread_args) {
            .n = n, .legal_states = legal_states, .legal_dict = legal_dict, .colorings = colorings,
            .start_idx = from, .end_idx = to, .result = result, .stop = stop, stop_after_first = stop_after_first, .num_done = num_done
        };

        if (force_threaded || num_threads > 1) {
            if (pthread_create(&pids[i],0,&orbit_thread_enter,&args[i])) {
                fprintf(stderr, "error: thread couldn't be created\n");
                exit(1);
            }
        } else { // single-threaded: don't create a new thread but perform the calculation here
            orbit_thread_enter(args);
        }
    }

    // Create calculation to track the progress
    legal_orbits_calculation calc = {
        .n = n, .num_threads = num_threads, .pids = pids, .args = args, .stop = stop, .num_colorings = colorings.len, .begin_ms = millis()
    };

    return calc;
}

// Find all legal orbits from the given legal states in the given coloring.
// The returned list contains a single legal state per legal orbit.
// When stop_after_first is true, set stop to true when a legal orbit was found.
static inline legal_orbits_result find_legal_orbits_single(int n, int* coloring, arr2d_fixed legal_states, bool* legal, legal_orbits_result append_result, bool* stop, bool stop_after_first) {
    legal_orbits_result result = append_result;
    bool found_orbits = false;

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

        if (orbit_legal) {
            if (found_orbits)
                result.states = append_arrv_single_into_last_row(result.states, state);
            else
                result.states = append_arrv_single(result.states, state);
            found_orbits = true;

            if (stop_after_first) {
                *stop = true; // signal other threads
                break;
            }
        }
    }

    if (found_orbits)
        result.colorings = append_arrf(result.colorings, coloring);

    return result;
}


/*
// Alternative algorithm: reduce all states into a single coset, then sort and analyze this list.
// We need a very quick (radix) sort for this to be efficient!
static inline arr2d_fixed find_legal_orbits_single_alternative(int n, int* coloring, arr2d_fixed legal_states, arr2d_fixed append_result) {
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