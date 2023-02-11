///
/// Unifies the qsort_r and qsort_s functions which are available on different systems.
/// Taken and modified from https://github.com/noporpoise/sort_r
///

/* Isaac Turner 29 April 2014 Public Domain */
#ifndef SORT_R_H_
#define SORT_R_H_

#include <stdlib.h> /* qsort_r(), qsort_s() */

/*

sort_r function to be exported.

Parameters:
  base is the array to be sorted
  nel is the number of elements in the array
  width is the size in bytes of each element of the array
  compar is the comparison function
  arg is a pointer to be passed to the comparison function

void sort_r(void *base, size_t nel, size_t width,
            int (*compar)(const void *_a, const void *_b, void *_arg),
            void *arg);

*/

#define _SORT_R_INLINE inline

#if (defined __APPLE__ || defined __MACH__ || defined __DARWIN__ || \
     (defined __FreeBSD__ && !defined(qsort_r)) || defined __DragonFly__)
#  define _SORT_R_BSD
#elif (defined __GLIBC__ || (defined (__FreeBSD__) && defined(qsort_r)))
#  define _SORT_R_LINUX
#elif (defined _WIN32 || defined _WIN64 || defined __WINDOWS__ || \
       defined __MINGW32__ || defined __MINGW64__)
#  define _SORT_R_WINDOWS
#  undef _SORT_R_INLINE
#  define _SORT_R_INLINE __inline
#endif

#if defined _SORT_R_BSD

/* Ensure qsort_r is defined */
extern void qsort_r(void *base, size_t nel, size_t width, void *thunk,
                    int (*compar)(void *_thunk,
                                    const void *_a, const void *_b));

#endif

#if defined _SORT_R_BSD || defined _SORT_R_WINDOWS

/* BSD (qsort_r), Windows (qsort_s) require argument swap */

struct sort_r_data
{
    void *arg;
    int (*compar)(const void *_a, const void *_b, void *_arg);
};

static _SORT_R_INLINE int sort_r_arg_swap(void *s,
                                            const void *a, const void *b)
{
    struct sort_r_data *ss = (struct sort_r_data*)s;
    return (ss->compar)(a, b, ss->arg);
}

#endif

#if defined _SORT_R_LINUX

typedef int(* __compar_d_fn_t)(const void *, const void *, void *);
extern void (qsort_r)(void *base, size_t nel, size_t width,
                        __compar_d_fn_t __compar, void *arg)
    __attribute__((nonnull (1, 4)));

#endif

/* implementation */

static _SORT_R_INLINE void sort_r(void *base, size_t nel, size_t width,
                                int (*compar)(const void *_a,
                                                const void *_b, void *_arg),
                                void *arg) {
#if defined _SORT_R_LINUX

    #if defined __GLIBC__ && ((__GLIBC__ < 2) || (__GLIBC__ == 2 && __GLIBC_MINOR__ < 8))

    fprintf(stderr, "there seems to be no qsort_s function on your system\n");
    exit(1);

    #else

    qsort_r(base, nel, width, compar, arg);

    #endif

#elif defined _SORT_R_BSD

    struct sort_r_data tmp;
    tmp.arg = arg;
    tmp.compar = compar;
    qsort_r(base, nel, width, &tmp, sort_r_arg_swap);

#elif defined _SORT_R_WINDOWS

    struct sort_r_data tmp;
    tmp.arg = arg;
    tmp.compar = compar;
    qsort_s(base, nel, width, sort_r_arg_swap, &tmp);

#else

    fprintf(stderr, "there seems to be no qsort_s function on your system\n");
    exit(1);

#endif
}

#undef _SORT_R_INLINE
#undef _SORT_R_WINDOWS
#undef _SORT_R_LINUX
#undef _SORT_R_BSD

#endif /* SORT_R_H_ */