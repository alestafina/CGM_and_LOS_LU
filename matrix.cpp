#include "laba3.h"

void matrix::input_matrix() {
   ifstream kuslau_,ig_, jg_, ggl_, ggu_, di_, pr_;
   kuslau_.open("BM/kuslau.txt");
   kuslau_ >> N >> eps >> maxiter;
   kuslau_.close();

   ig_.open("BM/ig.txt");
   ig.resize(N + 1);
   for (int i = 0; i < N + 1; ++i) {
      ig_ >> ig[i];
      ig[i]--;
   }
   ig_.close();

   jg_.open("BM/jg.txt");
   jg.resize(ig.back());
   for (int i = 0; i < jg.size(); ++i) {
      jg_ >> jg[i];
      jg[i]--;
   }
   jg_.close();

   ggl_.open("BM/ggl.txt");
   ggl.resize(ig.back());
   for (int i = 0; i < ggl.size(); ++i) {
      ggl_ >> ggl[i];
   }
   ggl_.close();
   
   ggu_.open("BM/ggu.txt");
   ggu.resize(ig.back());
   for (int i = 0; i < ggu.size(); ++i) {
      ggu_ >> ggu[i];
   }
   ggu_.close();
   
   di_.open("BM/di.txt");
   di.resize(N);
   for (int i = 0; i < N; ++i) {
      di_ >> di[i];
   }
   di_.close();
   
   pr_.open("BM/pr.txt");
   pr.resize(N);
   for (int i = 0; i < N; i++) {
      pr_ >> pr[i];
   }
   pr_.close();
}

void matrix::set_vector() {
   x.resize(N);
   x[0] = 1.0;
}

void matrix::output() {
   ofstream out("answer.txt");
   out.precision(15);
   for (int i = 0; i < N; i++)
      out << x[i] << endl;
   out.close();
}

vector<double> matrix::matr_vec_mult(vector<double> &x, bool flag) {
   vector<double> result(x.size()), L = ggl, U = ggu;
   if (flag) swap(L, U);
   for (int i = 0; i < x.size(); ++i) {
      result[i] = di[i] * x[i];
      for (int j = ig[i]; j < ig[i + 1]; j++) {
         result[i] += L[j] * x[jg[j]];
         result[jg[j]] += U[j] * x[i];
      }
   }
   return result;
}

double matrix::dot_product(vector<double> &a, vector<double> &b) {
   double result = 0.0;
   for (int i = 0; i < N; ++i)
      result += a[i] * b[i];
   return result;
}

void matrix::CGM() {
   double residual = 10.0, alpha, beta, r_prev, r_next, norm_b;
   set_vector();
   norm_b = dot_product(pr, pr);
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; i++) {
      r[i] = pr[i] - buf1[i];
   }
   r = matr_vec_mult(r, 1); // получили r0 для несимметричной матрицы
   z = r;
   while (k < maxiter && residual > eps) {
      buf1 = matr_vec_mult(z, 0);
      buf1 = matr_vec_mult(buf1, 1); // A'Az0
      alpha = dot_product(r, r) / dot_product(buf1, z);
      for (int i = 0; i < N; ++i)
         x[i] += alpha * z[i];
      r_prev = dot_product(r, r); // r(k-1)
      for (int i = 0; i < N; ++i) {
         r[i] -= alpha * buf1[i];
      }
      r_next = dot_product(r, r); // r(k)
      beta = r_next / r_prev;
      for (int i = 0; i < N; ++i)
         z[i] = r[i] + beta * z[i];
      residual = sqrt(r_next / norm_b);
      ++k;
   }
   cout << "CGM without precondition worked in " << k << " iterations. Residual is " << residual << endl;
}

void matrix::CGM_precond_diag() {
   double residual = 10.0, alpha, beta, r_prev, r_next, norm_b;
   vector<double> M (N); // M = diag(A'A) - матрица предобуславливания
   set_vector();
   for (int i = 0; i < N; ++i) {
      double sum = 0.0;
      for (int j = ig[i]; j < ig[i + 1]; j++) 
         sum += ggu[j] * ggu[j];
      for (int j = i + 1; j < N; j++)
         for (int q = ig[j]; q < ig[j + 1]; q++)
            if (jg[q] == i) sum += ggl[q] * ggl[q];
      M[i] = sum + di[i] * di[i];
   }

   norm_b = dot_product(pr, pr);
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; i++) {
      r[i] = pr[i] - buf1[i];
   }
   r = matr_vec_mult(r, 1); // получили r0 для несимметричной матрицы
   z.resize(N);
   for (int i = 0; i < N; i++) {
      z[i] = r[i] / M[i]; // z = M^(-1) * r0
   }
   while (k < maxiter && residual > eps) {
      buf2 = matr_vec_mult(z, 0);
      buf2 = matr_vec_mult(buf2, 1); // A'Az
      for (int i = 0; i < N; ++i)
         buf1[i] = r[i] / M[i]; // M^(-1) * r(k-1)
      alpha = dot_product(buf1, r) / dot_product(buf2, z);
      r_prev = dot_product(buf1, r); // r(k-1)
      for (int i = 0; i < N; ++i)
         x[i] += alpha * z[i];
      for (int i = 0; i < N; ++i) {
         r[i] -= alpha * buf2[i];
      }
      for (int i = 0; i < N; ++i)
         buf1[i] = r[i] / M[i]; // M^(-1) * r(k)
      r_next = dot_product(buf1, r); // r(k)
      beta = r_next / r_prev;
      for (int i = 0; i < N; ++i)
         z[i] = buf1[i] + beta * z[i];
      residual = sqrt(r_next / norm_b);
      ++k;
   }
   cout << "CGM with diagonal precondition worked in " << k << " iterations. Residual is " << residual << endl;
   
}

void matrix::LOS() {
   double residual = 10.0, alpha, beta, pp;
   set_vector();
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; i++) {
      r[i] = pr[i] - buf1[i];
   }
   residual = dot_product(r, r);
   z = r;
   p = matr_vec_mult(z, 0);
   while (k < maxiter && residual > eps) {
      pp = dot_product(p, p);
      alpha = dot_product(p, r) / pp;
      for (int i = 0; i < N; ++i) {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      buf1 = matr_vec_mult(r, 0);
      beta = -(dot_product(p, buf1) / pp);
      for (int i = 0; i < N; ++i) {
         z[i] = r[i] + beta * z[i];
         p[i] = buf1[i] + beta * p[i];
      }
      residual -= alpha * alpha * pp;
      ++k;
   }
   cout << "LOS without precondition worked in " << k << " iterations. Residual is " << residual << endl;
}

void matrix::LOS_precond_diag() { 
   // М = diag(LU) = diag(A) - матрица диагонального прдобуславливания,
   // diag(L) - содержит диагональные элементы А, 
   // diag(U) - содержит единицы на диагонали
   double residual = 10.0, alpha, beta, pp;
   set_vector();
   buf1 = matr_vec_mult(x, 0);
   r.resize(N);
   for (int i = 0; i < N; ++i) { // r0 = M(-1) * (f - A * x0)
      r[i] = (pr[i] - buf1[i]) / di[i];
   }
   residual = dot_product(r, r);
   z = r; // так как матрица U единичная
   p = matr_vec_mult(z, 0);
   for (int i = 0; i < N; ++i) { // p0 = M(-1) * A * z0
      p[i] /= di[i];
   }
   while (k < maxiter && residual > eps) {
      pp = dot_product(p, p);
      alpha = dot_product(p, r) / pp;
      for (int i = 0; i < N; ++i) {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      buf1 = matr_vec_mult(r, 0); // buf1 = A * r(k)
      for (int i = 0; i < N; ++i) {  // buf1 =  M(-1) * buf1
         buf1[i] /= di[i];
      }
      beta = -(dot_product(p, buf1) / pp);
      for (int i = 0; i < N; ++i) {
         z[i] = r[i] + beta * z[i];
         p[i] = buf1[i] + beta * p[i];
      }
      residual -= alpha * alpha * pp;
      ++k;
   }
   cout << "LOS with diagonal precondition worked in " << k << " iterations. Residual is " << residual << endl;
}

