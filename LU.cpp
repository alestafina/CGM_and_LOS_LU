#include "laba3.h"

void LU::ILU_precond() {
   ggl_LU = ggl;
   ggu_LU = ggu;
   di_LU = di;
   double sumL = 0.0, sumU = 0.0, sumD = 0.0;
   for (int k = 0; k < N; k++) {
      sumD = 0.0;
      int i0 = ig[k], i1 = ig[k + 1];
      int i = i0;
      for (; i0 < i1; i0++) {
         sumL = 0.0;
         sumU = 0.0;
         int j0 = i, j1 = i0;
         for (; j0 < j1; j0++) {
            int t0 = ig[jg[i0]], t1 = ig[jg[i0] + 1];
            for (; t0 < t1; t0++) {
               if (jg[j0] == jg[t0]) {
                  sumL += ggl_LU[j0] * ggu_LU[t0];
                  sumU += ggu_LU[j0] * ggl_LU[t0];
               }
            }
         }
         ggl_LU[i0] -= sumL;
         ggu_LU[i0] -= sumU;
         ggu_LU[i0] /= di_LU[jg[i0]];
         sumD += ggl_LU[i0] * ggu_LU[i0];
      }
      di_LU[k] -= sumD;
   }
}

vector<double> LU::direct(vector<double> &L, vector<double> &b) { // обрезанный прямой U(T)x = b
   vector<double> result(N);
   for (int i = 0; i < N; i++) {
      double sum = 0.0;
      for (int k0 = ig[i]; k0 < ig[i + 1]; k0++) {
         sum += result[jg[k0]] * L[k0];
      }
      result[i] = b[i] - sum;
   }
   return result;
}

vector<double> LU::direct(vector<double> &L, vector<double> &D, vector<double> &b) { // натуральный прямой Lx = b
   vector<double> result(N);
   for (int i = 0; i < N; i++) {
      double sum = 0.0;
      for (int k0 = ig[i]; k0 < ig[i + 1]; k0++) {
         sum += result[jg[k0]] * L[k0];
      }
      result[i] = (b[i] - sum) / D[i];
   }
   return result;
}

vector<double> LU::reverse(vector<double> &U, vector<double> &b) { // натуральный обратный Ux = b
   vector<double> result = b;
   for (int i = N - 1; i >= 0; i--) {
      int k0 = ig[i], k1 = ig[i + 1];
      int j;
      for (; k0 < k1; k0++) {
         j = jg[k0];
         result[j] -= result[i] * U[k0];
      }
   }
   return result;
}

vector<double> LU::reverse(vector<double> &U, vector<double> &D, vector<double> &b) { // обрезанный обратный L(T)x = b
   vector<double> result = b;
   for (int i = N - 1; i >= 0; i--) {
      int k0 = ig[i], k1 = ig[i + 1];
      int j;
      result[i] /= D[i];
      for (; k0 < k1; k0++) {
         j = jg[k0];
         result[j] -= result[i] * U[k0];
      }
   }
   return result;
}

void LU::CGM_precond_ILU() {
   double residual = 10.0, alpha, beta, norm_b, r_next, r_prev;
   set_vector();
   ILU_precond();
   norm_b = dot_product(pr, pr);
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; i++) { // r0 = f - Ax0
      r[i] = pr[i] - buf1[i];
   }
   r = direct(ggl_LU, di_LU, r); // r0 = L(-1)(f - Ax0)
   r = reverse(ggl_LU, di_LU, r); // r0 = L(-T)L(-1)(f - Ax0)
   r = matr_vec_mult(r, 1); // r0 = A(T)L(-T)L(-1)(f - Ax0)
   r = direct(ggu_LU, r); // r0 = U(-T)A(T)L(-T)L(-1)(f - Ax0)
   z = r;
   buf1.clear();
   buf1.resize(N);
   for (int i = 0; i < N; i++) { // x = Ux
      buf1[i] += x[i];
      for (int j = ig[i]; j < ig[i + 1]; j++)
         buf1[jg[j]] += ggu_LU[j] * x[i];			
   }
   x = buf1;
   while (k < maxiter && residual > eps) {
      r_prev = dot_product(r, r); // r(k-1)
      buf2 = reverse(ggu_LU, z);  // U(-1)z
      buf2 = matr_vec_mult(buf2, 0); // AU(-1)z
      buf2 = direct(ggl_LU, di_LU, buf2); // L(-1)AU(-1)z
      buf2 = reverse(ggl_LU, di_LU, buf2); // L(-T)L(-1)AU(-1)z
      buf2 = matr_vec_mult(buf2, 1); // A(T)L(-T)L(-1)AU(-1)z
      buf2 = direct(ggu_LU, buf2); // U(-T)A(T)L(-T)L(-1)AU(-1)z
      alpha = r_prev / dot_product(buf2, z);
      for (int i = 0; i < N; ++i) {
         x[i] += alpha * z[i];
         r[i] -= alpha * buf2[i];
      }
      r_next = dot_product(r, r);
      beta = r_next / r_prev;
      for (int i = 0; i < N; ++i) {
         z[i] = r[i] + beta * z[i];
      }
      residual = sqrt(r_next / norm_b);
      ++k;
   }
   x = reverse(ggu_LU, x); // x = U(-1)x
   cout << "CGM with ILU precondition worked in " << k << " iterations. Residual is " << residual << endl;
}



void LU::LOS_precond_ILU() {
   double residual = 10.0, alpha, beta, pp;
   set_vector();
   ILU_precond();
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; i++) {
      r[i] = pr[i] - buf1[i];
   }
   r = direct(ggl_LU, di_LU, r);
   z = reverse(ggu_LU, r);
   buf1 = matr_vec_mult(z, 0);  // Az
   p = direct(ggl_LU, di_LU, buf1); // L(-1)Az
   residual = dot_product(r, r);
   while (k < maxiter && residual > eps) {
      pp = dot_product(p, p);  
      alpha = dot_product(p, r) / pp;
      for (int i = 0; i < N; ++i) {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      buf2 = reverse(ggu_LU, r);  // U(-1)r
      buf1 = matr_vec_mult(buf2, 0); // AU(-1)r
      buf1 = direct(ggl_LU, di_LU, buf1); // L(-1)AU(-1)r
      beta = -(dot_product(p, buf1) / pp);
      for (int i = 0; i < N; ++i) {
         z[i] = buf2[i] + beta * z[i];
         p[i] = buf1[i] + beta * p[i];
      }
      residual -= alpha * alpha * pp;
      k++;
   }
   cout << "LOS with ILU precondition worked in " << k << " iterations. Residual is " << residual << endl;
}

