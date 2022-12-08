#include "laba3.h"

void Hilbert::set_Hilbert(int size) {
   N = size;
   int num = N * (N - 1) / 2;
   di.resize(size);
   for (int i = 0; i < size; i++) {
      di[i] = 1.0 / (2 * i + 1);
   }
   ig.resize(size + 1);
   ggl.resize(num);
   ggu.resize(num);
   jg.resize(num);
   ig[0] = 0;
   for (int i = 0; i < size; i++) {
      ig[i + 1] = ig[i] + i;
      for (int j = 0; j < i; j++) {
         ggl[ig[i] + j] = 1.0 / (i + j + 1);
         ggu[ig[i] + j] = ggl[ig[i] + j];
         jg[ig[i] + j] = j;
      }
   }
   set_pr();
}

void Hilbert::set_pr() {
   eps = 1e-15;
   maxiter = 1000;
   vector<double> vec(N);
   for (int i = 0; i < N; i++)
      vec[i] = i + 1;
   pr = matr_vec_mult(vec, 0);
}

