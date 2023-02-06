#ifndef UTILS_H
#define UTILS_H

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#include <stdio.h>
#include <stdbool.h>

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
    int len_guess; // a guess for total len that is used when successively filling the array without knowing the total rows
} arr2d_fixed;
arr2d_fixed arr2d_fixed_create_empty(int row_len, int len_guess);
arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len);
int get_arrf(arr2d_fixed arr, int i, int j); // first index by row, then by column
int get_arrf1d(arr2d_fixed arr, int i); // get the i'th entry of the contiguous data
arr2d_fixed append_arrf_single(arr2d_fixed arr, int val); // append a single value, given that row_len=1
arr2d_fixed append_arrf(arr2d_fixed arr, int* src); // reallocs the data to a larger size if len_guess is exceeded
void free_arrf(arr2d_fixed arr);
void print_arrf(arr2d_fixed arr);

// An arr2d_var is a 2D array of ints whose rows have not necessarily the same length.
typedef struct {
    int* data;
    int* start_indices; // has len+1 entries, the first entry being 0
    int len;
} arr2d_var;
arr2d_var arr2d_var_create_empty(int max_total_elems, int max_len);
arr2d_var arr2d_var_create_from(int* data, int* start_indices, int len);
int size_arrv(arr2d_var arr, int i);
int get_arrv(arr2d_var arr, int i, int j); // first index by row, then by column
arr2d_var append_arrv_single(arr2d_var arr, int val); // append a single value
arr2d_var append_arrv(arr2d_var arr, int* src, int n);
void free_arrv(arr2d_var arr);
void print_arrv(arr2d_var arr);

/******** COMBINATORICS ********/

int choose(int n, int k);
int ordered_choose(int n, int k);
arr2d_fixed do_choose(int n, int k, int* ptr);
arr2d_fixed do_ordered_choose(int n, int k, int* ptr);

#endif