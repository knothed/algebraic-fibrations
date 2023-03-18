///
/// Code that is used to check a stream of graphs (from geng) for fibering.
///

#include <stdio.h>
#include <stdbool.h>
#include <pthread.h>
#include "functionality.h"
#include "utils.h"

// Convert from geng's graph6 format to an adjacency matrix.
void read_adj_matrix_graph6(char* geng, int* adj) {
    int n = (int)(geng[0]-63);
    memset(adj,0,n*n*sizeof(int));

    char idx1 = 1;
    char idx2 = 0;
    char chr = 0;
    chr = geng[idx1] - 63;
    for (int i=1; i<n; i++) {
        for (int j=0; j<i; j++) {
            int v = (int)((chr >> (5-idx2)) & 1);
            adj[n*i+j] = v;
            adj[n*j+i] = v;

            idx2++;
            if (idx2 == 6) {
                idx2 = 0; idx1++;
                chr = geng[idx1] - 63;
            }
        }
    }
}

char* graph6_from_adj_matrix(arr2d_fixed adj) {
    int n = adj.len;
    char* res = malloc(3 + n*(n-1)/12);
    res[0] = n+63;

    char curr = 0;
    char idx1 = 1;
    char idx2 = 0;
    for (int i=1; i<n; i++) {
        for (int j=0; j<i; j++) {
            curr += get_arrf(adj,i,j) << (5-idx2);

            idx2++;
            if (idx2 == 6) {
                res[idx1] = curr + 63;
                idx2 = 0; idx1++;
                curr = 0;
            }
        }
    }

    if (idx2 > 0) {
        res[idx1] = curr + 63;
        idx1++;
    }
    res[idx1] = '\0';

    return res;
}

// A queue on which a stream of fibering calculations is performed.
typedef struct {
    pthread_t pid;
    int capacity;
    int start; // set by the worker thread
    int end; // set by the scheduler
    bool stop; // set by the scheduler
    int threads;

    arr2d_fixed* adj_data;
    arr2d_var* cliques_data;

    int checked_count;
    arr2d_fixed results; // all fibering graphs, consecutively

    FILE* results_file;
    pthread_mutex_t* mutex;

    uint64_t work_time;
    uint64_t wait_time;
} fibering_queue;

typedef struct {
    int n;
    int num_queues;
    fibering_queue** queues;
    uint64_t* wait_time; // all queues are full
    int64_t creation_time;
} fibering_scheduler;

// Queue loop running on a thread: wait for graphs coming in the queue and process them
void* queue_run(void* arg) {
    fibering_queue* queue = (fibering_queue*)arg;

    while (true) {
        if (queue->start == queue->end) {
            // wait a bit.
            // not too short, else waiting wastes CPU cycles
            // not too long, else waiting time adversely affects total speed
            delay(3);
            queue->wait_time += 3;

            // exit the loop only when all remaining tasks have been processed
            if (queue->stop && queue->start == queue->end)
                break;
        }
        else {
            int64_t t = millis();
            arr2d_fixed adj = queue->adj_data[queue->start];
            if (true) {
                queue->checked_count++;
                legal_orbits_result orbits = graph_fiberings(adj, queue->cliques_data[queue->start], 0, 0, false, false, queue->threads, true);
                bool fibers = orbits.colorings.len > 0;
                free_arrf(orbits.colorings);
                free_arrv(orbits.states);

                if (fibers) {
                    queue->results = append_arrf_multiple(queue->results, adj);
                    if (queue->results_file > 0) {
                        char* str = graph6_from_adj_matrix(adj);
                        pthread_mutex_lock(queue->mutex);
                        fputs(str, queue->results_file);
                        fputs("\n", queue->results_file);
                        fflush(queue->results_file);
                        pthread_mutex_unlock(queue->mutex);
                    }
                }
            }

            free_arrf(adj);
            free_arrv(queue->cliques_data[queue->start]);
            queue->work_time += millis() - t;

            // move on to the next graph
            queue->start++;
            if (queue->start == queue->capacity)
                queue->start = 0;
        }
    }
}

// Create num_threads worker threads waiting for graphs to process. All graphs must have the same number of vertices.
fibering_scheduler make_scheduler(int n, int num_queues, int capacity_per_queue, int threads_per_queue, char* results_file_str) {
    fibering_scheduler result;
    result.n = n;
    result.queues = malloc(num_queues*sizeof(fibering_queue*));
    result.num_queues = num_queues;
    result.wait_time = malloc(sizeof(uint64_t));
    *result.wait_time = 0;
    result.creation_time = millis();
    int capacity = capacity_per_queue + 1; // implementation detail

    // Create file & mutex
    pthread_mutex_t* mutex = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mutex, NULL);
    FILE* results_file = 0;
    if (strlen(results_file_str)) {
        results_file = fopen(results_file_str, "a");
        if (results_file <= 0) {
            fprintf(stderr, "File %s couldn't be opened\n", results_file_str);
            exit(1);
        }
    }

    for (int i=0; i<num_queues; i++) {
        fibering_queue* queue = malloc(sizeof(fibering_queue));
        queue->capacity = capacity;
        queue->threads = threads_per_queue;
        queue->start = 0;
        queue->end = 0;
        queue->stop = false;
        queue->adj_data = malloc(capacity*sizeof(arr2d_fixed));
        queue->cliques_data = malloc(capacity*sizeof(arr2d_var));
        queue->checked_count = 0;
        queue->results = arr2d_fixed_create_empty(n,n);
        queue->results_file = results_file;
        queue->mutex = mutex;
        queue->work_time = 0;
        queue->wait_time = 0;

        if (pthread_create(&queue->pid,0,&queue_run,queue)) {
            fprintf(stderr, "error: thread couldn't be created\n");
            exit(1);
        }

        result.queues[i] = queue;
    }

    return result;
}

// Add the graph to a free queue in the scheduler.
void add_to_scheduler(fibering_scheduler scheduler, arr2d_fixed adj, arr2d_var cliques) {
    // Find free queue
    int queue_idx = -1;
    while (queue_idx < 0) {
        for (int i=0; i<scheduler.num_queues; i++) {
            fibering_queue* queue = scheduler.queues[i];
            if ((queue->end+1)%queue->capacity != queue->start) {
                queue_idx = i;
                goto out;
            }
        }
        delay(3);
        (*scheduler.wait_time) += 3;
    }
    out: if (0) {}

    // Put work onto queue
    fibering_queue* queue = scheduler.queues[queue_idx];
    queue->adj_data[queue->end] = adj;
    queue->cliques_data[queue->end] = cliques;

    // advance queue_end
    queue->end++;
    if (queue->end == queue->capacity)
        queue->end = 0;
    if (queue->end == queue->start) {
        printf("Shouldn't happen!\n");
    }
}

typedef struct {
    int num_checked; // graphs checked for fibering
    int num_fiber;
    arr2d_fixed results;
} stream_result;

// Signal the scheduler that no new graphs will come in.
// Wait for all the worker threads to finish; then print work stats and return the result.
stream_result scheduler_finish(fibering_scheduler scheduler) {
    // Wait for queues to finish
    for (int i=0; i<scheduler.num_queues; i++)
        scheduler.queues[i]->stop = true;
    for (int i=0; i<scheduler.num_queues; i++)
        pthread_join(scheduler.queues[i]->pid, 0);

    // Collect result
    stream_result result;
    result.num_checked = 0;
    result.results = arr2d_fixed_create_empty(scheduler.n, scheduler.n);
    uint64_t total_work_time = 0;

    for (int i=0; i<scheduler.num_queues; i++) {
        result.num_checked += scheduler.queues[i]->checked_count;
        result.results = append_arrf_multiple(result.results, scheduler.queues[i]->results);
        total_work_time += scheduler.queues[i]->work_time;
    }
    result.num_fiber = result.results.len / scheduler.n;

    if (scheduler.queues[0]->results_file > 0)
        fclose(scheduler.queues[0]->results_file);

    // Print stats
    printf("\nStreaming finished. Graphs checked: %d, %d of which fiber(s).\n", result.num_checked, result.num_fiber);
    char* total = pretty_ms(millis() - scheduler.creation_time, true);
    printf("Took %s in total.\n", total);

    char* work = pretty_ms(total_work_time, true);
    printf(" • raw search time: %s (queue distribution: ", work);
    for (int i=0; i<scheduler.num_queues; i++) {
        char* time = pretty_ms(scheduler.queues[i]->work_time, true);
        printf("%s", time);
        if (i<scheduler.num_queues-1) printf(", ");
        free(time);
    }
    printf(")\n");

    char* full = pretty_ms(*scheduler.wait_time, true);
    printf(" • all queues full: %s\n", full);
    printf(" • queues empty: ");
    for (int i=0; i<scheduler.num_queues; i++) {
        char* time = pretty_ms(scheduler.queues[i]->wait_time, true);
        printf("%s", time);
        if (i<scheduler.num_queues-1) printf(", ");
        free(time);
    }
    printf("\n");

    free(total);
    free(work);
    free(full);

    // Free stuff
    free(scheduler.queues[0]->mutex);
    for (int i=0; i<scheduler.num_queues; i++) {
        free(scheduler.queues[i]->adj_data);
        free(scheduler.queues[i]->cliques_data);
        free(scheduler.queues[i]);
    }
    free(scheduler.queues);
    free(scheduler.wait_time);

    return result;
}