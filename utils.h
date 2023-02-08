#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <time.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static inline int log2_int(int a) {
    if (a <= 0) return -1;
    int r = 0;
    int b = a;
    while (b >>= 1) r++;
    return r;
}

int64_t millis();

/******** 2D ARRAY OF INTS ********/
// We define two types of 2D arrays, i.e. arrays of arrays of ints:
// arr2d_fixed is a 2D array where all rows have a fixed, predefined length.
// arr2d_var is a 2D array where the rows can have different arbitrary lengths.
// In both cases, the number of columns is variable and may change.
// Both are just convenient wrappers around a contiguous memory segment.

// An arr2d_fixed is a 2D array of ints whose rows all have the same length.
typedef struct {
    int* data; // in total row_len * len entries
    int row_len;
    int len;
    int capacity; // a guess for total len that is used when successively filling the array without knowing the total rows
} arr2d_fixed;
arr2d_fixed arr2d_fixed_create_empty(int row_len, int capacity);
arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len);
arr2d_fixed append_arrf_multiple(arr2d_fixed arr, arr2d_fixed new); // append multiple rows to the array
void free_arrf(arr2d_fixed arr);
void print_arrf(arr2d_fixed arr);
void print_arrf_row(arr2d_fixed arr, int i);

// first index by row, then by column
static inline int get_arrf(arr2d_fixed arr, int i, int j) {
    return arr.data[i*arr.row_len+j];
}

// get the i'th entry of the contiguous data
static inline int get_arrf1d(arr2d_fixed arr, int i) {
    return arr.data[i];
}

// reallocs the data to a larger size if capacity is exceeded
static inline arr2d_fixed append_arrf(arr2d_fixed arr, int* src) {
    arr2d_fixed other = arr2d_fixed_create_from(src, arr.row_len, 1);
    return append_arrf_multiple(arr, other);
}

// append a single value, given that row_len=1
static inline arr2d_fixed append_arrf_single(arr2d_fixed arr, int val) {
    int p[1] = {val};
    return append_arrf(arr, p);
}

// An arr2d_var is a 2D array of ints whose rows have not necessarily the same length.
typedef struct {
    int* data;
    int* start_indices; // has len+1 entries, the first entry being 0
    int len;
} arr2d_var;
arr2d_var arr2d_var_create_empty(int max_total_elems, int max_len);
arr2d_var arr2d_var_create_from(int* data, int* start_indices, int len);
void free_arrv(arr2d_var arr);
void print_arrv(arr2d_var arr);
void print_arrv_row(arr2d_var arr, int i);

// Get the size of the i'th subarray.
static inline int size_arrv(arr2d_var arr, int i) {
    return arr.start_indices[i+1]-arr.start_indices[i];
}

// Get the j'th element of the i'th subarray.
static inline int get_arrv(arr2d_var arr, int i, int j) {
    return arr.data[arr.start_indices[i]+j];
}

static inline arr2d_var append_arrv(arr2d_var arr, int* src, int n) {
    arr2d_var new = arr;
    memcpy(arr.data+arr.start_indices[arr.len], src, n*sizeof(int));
    new.len++;
    new.start_indices[new.len] = new.start_indices[new.len-1]+n;
    return new;
}

// just append a single value
static inline arr2d_var append_arrv_single(arr2d_var arr, int val) {
    int p[1] = {val};
    return append_arrv(arr, p, 1);
}


/******** COMBINATORICS ********/

int choose(int n, int k);
int ordered_choose(int n, int k);
arr2d_fixed do_choose(int n, int k, int* ptr);
arr2d_fixed do_ordered_choose(int n, int k, int* ptr);

#endif