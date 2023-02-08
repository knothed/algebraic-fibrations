///
/// Utilities.
///

#include "utils.h"

int64_t millis() {
    struct timespec now;
    timespec_get(&now, TIME_UTC);
    return ((int64_t) now.tv_sec) * 1000 + ((int64_t) now.tv_nsec) / 1000000;
}

/******** ARR2D_FIXED ********/

arr2d_fixed arr2d_fixed_create_empty(int row_len, int capacity) {
    arr2d_fixed arr;
    arr.data = malloc(capacity*row_len*sizeof(int));
    arr.row_len = row_len;
    arr.len = 0;
    arr.capacity = capacity;
    return arr;
}

arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len) {
    arr2d_fixed arr;
    arr.data = data;
    arr.row_len = row_len;
    arr.len = len;
    arr.capacity = len;
    return arr;
}

// Assumes arr.row_len == other.row_len
arr2d_fixed append_arrf_multiple(arr2d_fixed arr, arr2d_fixed other) {
    arr2d_fixed new = arr;

    // realloc if the array is full
    new.len = arr.len + other.len;
    while (new.len > new.capacity)
        new.capacity = new.capacity + (new.capacity >> 1) + (new.capacity >> 3) + 1; // ca. phi
    if (new.capacity > arr.capacity)
        new.data = realloc(new.data, new.capacity*new.row_len*sizeof(int));

    memcpy(new.data+new.row_len*arr.len, other.data, other.len*other.row_len*sizeof(int));
    return new;
}

void free_arrf(arr2d_fixed arr) {
    free(arr.data);
}

void print_arrf(arr2d_fixed arr) {
    printf("{");
    for (int i=0; i<arr.len; i++) {
        if (i > 0) printf(", ");
        printf("{");
        for (int j=0; j<arr.row_len; j++) {
            if (j > 0) printf(",");
            printf("%d", get_arrf(arr,i,j));
        }
        printf("}");
    }
    printf("}\n");
}

void print_arrf_row(arr2d_fixed arr, int i) {
    printf("{");
    for (int j=0; j<arr.row_len; j++) {
        if (j > 0) printf(",");
        printf("%d", get_arrf(arr,i,j));
    }
    printf("}\n");
}

/******** ARR2D_VAR ********/

arr2d_var arr2d_var_create_empty(int max_total_elems, int max_len) {
    arr2d_var arr;
    arr.data = (int*)malloc(max_total_elems*sizeof(int));
    arr.start_indices = (int*)malloc((max_len+1)*sizeof(int));
    arr.len = 0;
    memset(arr.data,0,max_total_elems*sizeof(int));
    memset(arr.start_indices,0,(max_len+1)*sizeof(int));
    return arr;
}

arr2d_var arr2d_var_create_from(int* data, int* start_indices, int len) {
    arr2d_var arr;
    arr.data = data;
    arr.start_indices = start_indices;
    arr.len = len;
    return arr;
}

void free_arrv(arr2d_var arr) {
    free(arr.data);
    free(arr.start_indices);
}

void print_arrv(arr2d_var arr) {
    printf("{");
    for (int i=0; i<arr.len; i++) {
        if (i > 0) printf(", ");
        printf("{");
        for (int j=0; j<size_arrv(arr,i); j++) {
            if (j > 0) printf(",");
            printf("%d", get_arrv(arr,i,j));
        }
        printf("}");
    }
    printf("}\n");
}

void print_arrv_row(arr2d_var arr, int i) {
    printf("{");
    for (int j=0; j<size_arrv(arr,i); j++) {
        if (j > 0) printf(",");
        printf("%d", get_arrv(arr,i,j));
    }
    printf("}\n");
}


/******** COMBINATORICS ********/

void permute(int *subset, int start, int k, int* res, int* count);
void subset_helper(int *subset, int n, int k, int index, int start, int* res, int* count, bool ordered);
static inline void swap(int *a, int *b);

// Calculate (n choose k) * k!.
int ordered_choose(int n, int k) {
    int result = 1;
    for (int i=n; i>n-k; i--)
        result *= i;
    return result;
}

// Calculate (n choose k).
int choose(int n, int k) {
    if (k==0) return 1;
    return (n * choose(n-1, k-1)) / k;
}

// Fill the given array with all (n choose k) options of choosing k numbers of 0..<n in an unordered fashion.
// Return an arr2d_fixed which holds the given data pointer: for performance reasons, this array can be created on the caller stack.
arr2d_fixed do_choose(int n, int k, int* ptr) {
    int subset[k];
    int count = 0;
    subset_helper(subset,n,k,0,0,ptr,&count,false);
    return arr2d_fixed_create_from(ptr, k, choose(n,k));
}

// Fill the given array with all (n choose k) * k! options of choosing k numbers of 0..<n in an ordered fashion.
// Return an arr2d_fixed which holds the given data pointer: for performance reasons, this array can be created on the caller stack.
arr2d_fixed do_ordered_choose(int n, int k, int* ptr) {
    int subset[k];
    int count = 0;
    subset_helper(subset,n,k,0,0,ptr,&count,true);
    return arr2d_fixed_create_from(ptr, k, ordered_choose(n,k));
}

void permute(int *subset, int start, int k, int* res, int* count) {
    if (start == k) {
        memcpy(res+k*(*count),subset,k*sizeof(int));
        (*count)++;
        return;
    }

    for (int i=start; i<k; i++) {
        swap(subset+start, subset+i);
        permute(subset,start+1,k,res,count);
        swap(subset+start, subset+i);
    }
}

void subset_helper(int *subset, int n, int k, int index, int start, int* res, int* count, bool ordered) {
    if (index == k) {
        if (ordered) { // generate all permutations of subset
            permute(subset,0,k,res,count);
            return;
        } else { // add to result
            memcpy(res+k*(*count),subset,k*sizeof(int));
            (*count)++;
            return;
        }
    }

    for (int i=start; i<n; i++) {
        subset[index] = i;
        subset_helper(subset,n,k,index+1,i+1,res,count,ordered);
    }
}

static inline void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}