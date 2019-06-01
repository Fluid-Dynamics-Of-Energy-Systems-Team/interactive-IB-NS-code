


#ifndef DEFINES_H_
#define DEFINES_H_



#define MAX(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a >= _b ? _a : _b; })
#define MIN(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a <= _b ? _a : _b; })

#define ind(i,j)  ((i)*jmax+(j))
#define indp(i,j) ((i)*(jmax-2)+(j))

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }


#define NUM_THREADS 4
#define NUM_THREADS_FFT 4


#endif












