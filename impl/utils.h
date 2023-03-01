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

/******** THREAD DELAY ********/

#if defined(__WIN32__) || defined(_WIN32) || defined(WIN32) || defined(__WINDOWS__) || defined(__TOS_WIN__)

#include <windows.h>

inline void delay(unsigned long ms) {
    Sleep(ms);
}

#else  /* presume POSIX */

#include <unistd.h>

inline void delay(unsigned long ms) {
    usleep(ms * 1000);
}

#endif

/******** 2D ARRAY OF INTS ********/
// We define two types of 2D arrays, i.e. arrays of arrays of ints:
// arr2d_fixed is a 2D array where all rows have a fixed, predefined length.
// arr2d_var is a 2D array where the rows can have different arbitrary lengths.
// In both cases, the number of columns is variable and may change.
// Both are just convenient wrappers around a contiguous memory segment.

int phi_times(int x) {
    return x + (x >> 1) + (x >> 3) + 1;
}

// An arr2d_fixed is a 2D array of ints whose rows all have the same length.
typedef struct {
    int* data; // in total row_len * len entries
    int row_len;
    int len;
    int capacity; // size of data that is allocated
} arr2d_fixed;
arr2d_fixed arr2d_fixed_create_empty(int row_len, int capacity);
arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len);
arr2d_fixed append_arrf_multiple(arr2d_fixed arr, arr2d_fixed new); // append multiple rows to the array
#define free_arrf(...) free_arrfs(sizeof((arr2d_fixed[]) {__VA_ARGS__}) / sizeof(arr2d_fixed), __VA_ARGS__)
static void free_arrfs(int n, ...);
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
    int* end_indices; // has len entries, the last one being the total length
    int len;
    int num_rows_capacity; // size of end_indices that is allocated
    int total_capacity; // size of data that is allocated
} arr2d_var;
arr2d_var arr2d_var_create_empty(int total_capacity, int num_rows_capacity);
arr2d_var arr2d_var_create_from(int* data, int* end_indices, int len);
arr2d_var append_arrv_multiple(arr2d_var arr, arr2d_var other); // append multiple rows to the array
#define free_arrv(...) free_arrvs(sizeof((arr2d_var[]) {__VA_ARGS__}) / sizeof(arr2d_var), __VA_ARGS__)
static void free_arrvs(int n, ...);
void print_arrv(arr2d_var arr);
void print_arrv_row(arr2d_var arr, int i);

static inline int total_len_arrv(arr2d_var arr) {
    if (arr.len==0) return 0;
    return arr.end_indices[arr.len-1];
}

// Get the size of the i'th subarray.
static inline int size_arrv(arr2d_var arr, int i) {
    if (i==0) return arr.end_indices[0];
    return arr.end_indices[i]-arr.end_indices[i-1];
}

// Get the j'th element of the i'th subarray.
static inline int get_arrv(arr2d_var arr, int i, int j) {
    if (i==0) return arr.data[j];
    return arr.data[arr.end_indices[i-1]+j];
}

// reallocs the data to a larger size if capacity is exceeded
static inline arr2d_var append_arrv(arr2d_var arr, int* src, int n) {
    int index[1] = {n};
    arr2d_var other = arr2d_var_create_from(src, index, 1);
    return append_arrv_multiple(arr, other);
}

// append a new row with a single value
static inline arr2d_var append_arrv_single(arr2d_var arr, int val) {
    int p[1] = {val};
    return append_arrv(arr, p, 1);
}

// append a single value to the last row.
// precondition: len > 0
static inline arr2d_var append_arrv_single_into_last_row(arr2d_var arr, int val) {
    arr2d_var new = append_arrv_single(arr, val);
    new.len--;
    new.end_indices[new.len-1]++;
    return new;
}


/******** COMBINATORICS ********/

int choose(int n, int k);
int ordered_choose(int n, int k);
arr2d_fixed do_choose(int n, int k, int* ptr);
arr2d_fixed do_ordered_choose(int n, int k, int* ptr);


/******** PRETTY PRINTING ********/

// Convert a timespan, given in milliseconds, into a more human readable form.
char* pretty_ms(uint64_t ms, bool subsecond_precision);

// Equip a nonnegative integer with thousands delimiters.
char* pretty_int(int num);

// Print a progress bar onto the current line.
void print_progress(char* prefix, double progress, int64_t estimated_ms);

#endif