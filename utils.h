#ifndef UTILS_H
#define UTILS_H

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#include <stdio.h>
#include <stdbool.h>

/******** ARRAY OF ARRAYS OF DIFFERENT LENGTHS ********/

// An arr_of_arrs is a 2D array (of ints) whose entries have not necessarily the same length.
typedef struct {
    int* data;
    int* start_indices;
    int len;
} arr_of_arrs;

arr_of_arrs create_arr_of_arrs(int max_total_elems, int max_arrs);
arr_of_arrs comprise_arr_of_arrs(int* data, int* start_indices, int* len);

int arr_size(arr_of_arrs arr, int i);
int arr_get(arr_of_arrs arr, int i, int j);
arr_of_arrs arr_append(arr_of_arrs arr, int* src, int n);

void free_arr(arr_of_arrs arr);
void print_arr(arr_of_arrs arr);

/******** COMBINATORICS ********/

int choose(int n, int k);
int ordered_choose(int n, int k);
void do_choose(int n, int k, int* res);
void do_ordered_choose(int n, int k, int* res);

#endif