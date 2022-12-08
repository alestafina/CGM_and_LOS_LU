#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono> 

using namespace std;

class matrix {
   // N - ����������� �������, k - ����� ��������,
   // maxiter - ������������ ����� ��������,
   // x - ������� ������, r - ������ �������, 
   // z - ������ ������, p - ��������������� ������
protected:
   int N, maxiter, k;
   double eps;
   vector<double> ggl, ggu, di, pr, x;
   vector<double> buf1, buf2, r, z, p;
   vector<int> ig, jg;
public:
   matrix() {
      this->N = 0;
      this->maxiter = 0;
      this->k = 1;
      this->eps = 0.0;
      this->ggl = vector<double>();
      this->ggu = vector<double>();
      this->di = vector<double>();
      this->pr = vector<double>();
      this->x = vector<double>();
      this->r = vector<double>();
      this->z = vector<double>();
      this->p = vector<double>();
      this->ig = vector<int>();
      this->jg = vector<int>();
   }

   void input_matrix();
   void set_vector();
   void output();
   
   vector<double> matr_vec_mult(vector<double> &x, bool flag);
   double dot_product(vector<double> &a, vector<double> &b);

   void CGM();
   void CGM_precond_diag();
   void LOS();
   void LOS_precond_diag();
};

class LU : public matrix {
private:
   vector<double> ggl_LU, ggu_LU, di_LU;
public:
   void ILU_precond(); // �������� LU ������������

   vector<double> direct(vector<double> &L, vector<double> &b); // ��� ����������������� U
   vector<double> direct(vector<double> &L, vector<double> &D, vector<double> &b); // ��� ��������� L
   vector<double> reverse(vector<double> &U, vector<double> &b); // ��� ��������� U
   vector<double> reverse(vector<double> &U, vector<double> &D, vector<double> &b); // ��� ����������������� L
   
   void CGM_precond_ILU();
   void LOS_precond_ILU();
};

class Hilbert : public LU {
public:
   void set_Hilbert(int size);
   void set_pr();

};

