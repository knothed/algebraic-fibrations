///
/// Utilities.
///

#include "utils.h"
#include <stdarg.h>

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
        new.capacity = phi_times(new.capacity);
    if (new.capacity > arr.capacity)
        new.data = realloc(new.data, new.capacity*new.row_len*sizeof(int));

    memcpy(new.data+new.row_len*arr.len, other.data, other.len*other.row_len*sizeof(int));
    return new;
}

static void free_arrfs(int n, ...) {
    va_list args;
    va_start(args, n);
    while (n--) {
        arr2d_fixed arr = va_arg(args, arr2d_fixed);
        free(arr.data);
    }
    va_end(args);
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

arr2d_var arr2d_var_create_empty(int total_capacity, int num_rows_capacity) {
    arr2d_var arr;
    arr.data = (int*)malloc(total_capacity*sizeof(int));
    arr.end_indices = (int*)malloc(num_rows_capacity*sizeof(int));
    arr.len = 0;
    arr.total_capacity = total_capacity;
    arr.num_rows_capacity = num_rows_capacity;
    return arr;
}

arr2d_var arr2d_var_create_from(int* data, int* end_indices, int len) {
    arr2d_var arr;
    arr.data = data;
    arr.end_indices = end_indices;
    arr.len = len;
    arr.total_capacity = end_indices[len-1];
    arr.num_rows_capacity = len;
    return arr;
}

arr2d_var append_arrv_multiple(arr2d_var arr, arr2d_var other) {
    arr2d_var new = arr;
    new.len = arr.len + other.len;

    int total_arr_len = total_len_arrv(arr); // realloc of new.end_indices destroys this

    // realloc end_indices
    while (new.len > new.num_rows_capacity)
        new.num_rows_capacity = phi_times(new.num_rows_capacity);
    if (new.num_rows_capacity > arr.num_rows_capacity)
        new.end_indices = realloc(new.end_indices, new.num_rows_capacity*sizeof(int));

    for (int i=0; i<other.len; i++)
        new.end_indices[arr.len+i] = total_arr_len + other.end_indices[i];

    // realloc data
    while (total_len_arrv(new) > new.total_capacity)
        new.total_capacity = phi_times(new.total_capacity);
    if (new.total_capacity > arr.total_capacity)
        new.data = realloc(new.data, new.total_capacity*sizeof(int));

    memcpy(new.data+total_arr_len, other.data, total_len_arrv(other)*sizeof(int));
    return new;
}

static void free_arrvs(int n, ...) {
    va_list args;
    va_start(args, n);
    while (n--) {
        arr2d_var arr = va_arg(args, arr2d_var);
        free(arr.data);
        free(arr.end_indices);
    }
    va_end(args);
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

static inline void swap_char(char *a, char *b) {
    char tmp = *a;
    *a = *b;
    *b = tmp;
}


/******** PRETTY PRINTING ********/

// Convert a timespan, given in milliseconds, into a more human readable form.
char* pretty_ms(uint64_t ms, bool subsecond_precision) {
    char* result = malloc(8*sizeof(char));
    uint64_t s = (ms+500)/1000;

    if (subsecond_precision && ms+5 < 1000) {
        snprintf(result, 8, "%.2fs", ((double)ms+5)/1000.0);
    } else if (subsecond_precision && ms+50 < 10 * 1000) {
        snprintf(result, 8, "%.1fs", ((double)ms+50)/1000.0);
    } else if (s < 600) {
        snprintf(result, 8, "%llds", s);
    } else if (s < 600*600) {
        snprintf(result, 8, "%lldm", s/60);
    } else {
        snprintf(result, 8, "%lldh", s/3600);
    }

    return result;
}

// Equip a nonnegative integer with thousands delimiters.
char* pretty_int(int num) {
    int len = 3+log2_int(num)/2;
    char* result = malloc(len*sizeof(char));

    // Write string into buffer
    int i=0; int j=0;
    do {
        int digit = num%10;
        result[i+j] = digit+'0';
        num /= 10;

        i++;
        if (num > 0 && i%3==0) {
            result[i+j] = '\'';
            j++;
        }
    } while (num);

    // Reverse string
    len = i+j;
    for (int i=0; i<len/2; i++)
        swap_char(result+i, result+len-1-i);
    result[len] = 0;

    return result;
}