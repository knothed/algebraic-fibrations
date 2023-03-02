///
/// Code that interacts with geng to check a stream of graphs for fibering.
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
    for (int i=0; i<n; i++) {
        for (int j=i+1; j<n; j++) {
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

// A queue on which a stream of fibering calculations is performed.
typedef struct {
    pthread_t pid;
    int capacity;
    int start; // set by the worker thread
    int end; // set by the scheduler
    bool stop; // set by the scheduler

    arr2d_fixed* adj_data;
    arr2d_var* cliques_data;

    int hyp_count;
    int result_count; // set by the worker thread

    uint64_t work_time;
} fibering_queue;

typedef struct {
    int num_queues;
    fibering_queue** queues;
} fibering_scheduler;

// Queue loop running on a thread: wait for graphs coming in the queue and process them
void* queue_run(void* arg) {
    fibering_queue* queue = (fibering_queue*)arg;

    while (true) {
        if (queue->start == queue->end) {
            // wait a bit.
            // not too short, else waiting wastes CPU cycles
            // not too long, else waiting time adversely affects total speed
            delay(2);

            // exit the loop only when all remaining tasks have been processed
            if (queue->stop && queue->start == queue->end)
                break;
        }
        else {
            // printf("%d, %lld\n",queue->end-queue->start, millis()-t);
            int64_t t = millis();
            arr2d_fixed adj = queue->adj_data[queue->start];
            if (true) {
                queue->hyp_count++;
                legal_orbits_result orbits = graph_fiberings(adj, queue->cliques_data[queue->start], 0, 0, false, false, 3, true);
                bool fibers = orbits.colorings.len > 0;
                free_arrf(orbits.colorings);
                free_arrv(orbits.states);
                if (fibers)
                    queue->result_count++;
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

// Create num_threads worker threads waiting for graphs to process.
fibering_scheduler make_scheduler(int num_threads, int capacity_per_thread) {
    fibering_scheduler result;
    result.queues = malloc(num_threads*sizeof(fibering_queue*));
    result.num_queues = num_threads;

    for (int i=0; i<num_threads; i++) {
        fibering_queue* queue = malloc(sizeof(fibering_queue));
        queue->capacity = capacity_per_thread;
        queue->start = 0;
        queue->end = 0;
        queue->stop = false;
        queue->adj_data = malloc(capacity_per_thread*sizeof(arr2d_fixed));
        queue->cliques_data = malloc(capacity_per_thread*sizeof(arr2d_var));
        queue->hyp_count = 0;
        queue->result_count = 0;
        queue->work_time = 0;

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
        delay(5);
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

// Signal the scheduler that no new graphs will come in.
// Wait for all the worker threads to finish and return the result.
void scheduler_finish(fibering_scheduler scheduler) {
    // Wait for queues to finish
    for (int i=0; i<scheduler.num_queues; i++)
        scheduler.queues[i]->stop = true;
    for (int i=0; i<scheduler.num_queues; i++)
        pthread_join(scheduler.queues[i]->pid, 0);

    // Collect result
    int total_hyp = 0;
    int total_fiber = 0;
    uint64_t total_work_time = 0;
    for (int i=0; i<scheduler.num_queues; i++) {
        total_hyp += scheduler.queues[i]->hyp_count;
        total_fiber += scheduler.queues[i]->result_count;
        total_work_time += scheduler.queues[i]->work_time;
    }
    printf("finished. hyp: %d, fiber: %d, c-work: %lld\n", total_hyp, total_fiber, total_work_time);

    for (int i=0; i<scheduler.num_queues; i++) {
        printf("q %d: %lld ", i, scheduler.queues[i]->work_time);
        if (i<scheduler.num_queues-1) printf("; "); else printf("\n");
    }
}