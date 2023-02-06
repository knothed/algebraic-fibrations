///
/// Utilities.
///

#include "utils.h"

/******** ARR2D_FIXED ********/

arr2d_fixed arr2d_fixed_create_empty(int row_len, int len_guess) {
    arr2d_fixed arr;
    arr.data = malloc(len_guess*row_len*sizeof(int));
    arr.row_len = row_len;
    arr.len = 0;
    arr.len_guess = len_guess;
    return arr;
}

arr2d_fixed arr2d_fixed_create_from(int* data, int row_len, int len) {
    arr2d_fixed arr;
    arr.data = data;
    arr.row_len = row_len;
    arr.len = len;
    arr.len_guess = len;
    return arr;
}

int get_arrf(arr2d_fixed arr, int i, int j) {
    return arr.data[i*arr.row_len+j];
}

int get_arrf1d(arr2d_fixed arr, int i) {
    return arr.data[i];
}

arr2d_fixed append_arrf_single(arr2d_fixed arr, int val) {
    int p[1] = {val};
    return append_arrf(arr, p);
}

arr2d_fixed append_arrf(arr2d_fixed arr, int* src) {
    arr2d_fixed new = arr;

    // realloc if the array is full
    if (arr.len == arr.len_guess) {
        int new_size = arr.len_guess + (arr.len_guess >> 1) + (arr.len_guess >> 3) + 1; // ca. phi
        new.data = realloc(arr.data, new_size*new.row_len*sizeof(int));
        new.len_guess = new_size;
    }

    memcpy(new.data+new.row_len*new.len, src, new.row_len*sizeof(int));
    new.len++;
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

// Get the size of the i'th subarray.
int size_arrv(arr2d_var arr, int i) {
    return arr.start_indices[i+1]-arr.start_indices[i];
}

// Get the j'th element of the i'th subarray.
int get_arrv(arr2d_var arr, int i, int j) {
    return arr.data[arr.start_indices[i]+j];
}

arr2d_var append_arrv_single(arr2d_var arr, int val) {
    int p[1] = {val};
    return append_arrv(arr, p, 1);
}

arr2d_var append_arrv(arr2d_var arr, int* src, int n) {
    arr2d_var new = arr;
    memcpy(arr.data+arr.start_indices[arr.len], src, n*sizeof(int));
    new.len++;
    new.start_indices[new.len] = new.start_indices[new.len-1]+n;
    return new;
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


/******** COMBINATORICS ********/

void permute(int *subset, int start, int k, int* res, int* count);
void subset_helper(int *subset, int n, int k, int index, int start, int* res, int* count, bool ordered);
void swap(int *a, int *b);

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

void swap(int *a, int *b) {
    int tmp = *a;
    *a = *b;
    *b = tmp;
}