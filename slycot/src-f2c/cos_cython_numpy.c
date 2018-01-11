/*
2.8.5.2. Numpy Support, 2.8.5. Cython, http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id13
*/

#include <math.h>

/*  Compute the cosine of each element in in_array, storing the result in
 *  out_array. */
void cos_cython_numpy_func(double * in_array, double * out_array, int size){
    int i;
    for(i=0;i<size;i++){
        out_array[i] = cos(in_array[i]);
    }
}
