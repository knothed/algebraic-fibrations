///
/// Utilities.
///

#include "utils.h"

arr_of_arrs create_arr_of_arrs(int max_total_elems, int max_arrs) {
    arr_of_arrs arr;
    arr.data = (int*)malloc(max_total_elems*sizeof(int));
    arr.start_indices = (int*)malloc((max_arrs+1)*sizeof(int));
    arr.len = 0;
    memset(arr.data,0,max_total_elems*sizeof(int));
    memset(arr.start_indices,0,(max_arrs+1)*sizeof(int));
    return arr;
}

arr_of_arrs comprise_arr_of_arrs(int* data, int* start_indices, int* len) {
    arr_of_arrs arr;
    arr.data = data;
    arr.start_indices = start_indices;
    arr.len = len;
    return arr;
}

// Get the size of the i'th subarray.
int arr_size(arr_of_arrs arr, int i) {
    return arr.start_indices[i+1]-arr.start_indices[i];
}

// Get the j'th element of the i'th subarray.
int arr_get(arr_of_arrs arr, int i, int j) {
    return arr.data[arr.start_indices[i]+j];
}

arr_of_arrs arr_append(arr_of_arrs arr, int* src, int n) {
    arr_of_arrs new = arr;
    memcpy(arr.data+arr.start_indices[arr.len], src, n*sizeof(int));
    new.len++;
    new.start_indices[new.len] = new.start_indices[new.len-1]+n;
    return new;
}

void free_arr(arr_of_arrs arr) {
    free(arr.data);
    free(arr.start_indices);
}

void print_arr(arr_of_arrs arr) {
    printf("{");
    for (int i=0; i<arr.len; i++) {
        if (i > 0) printf(", ");
        printf("{");
        for (int j=0; j<arr_size(arr,i); j++) {
            if (j > 0) printf(",");
            printf("%d", arr_get(arr,i,j));
        }
        printf("}");
    }
    printf("}\n");
}