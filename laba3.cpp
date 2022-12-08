#include "laba3.h"

int main() {
   LU A;
   A.input_matrix();
   auto start = chrono::system_clock::now();
   A.CGM();
   auto end = chrono::system_clock::now();
   cout << chrono::duration_cast<chrono::microseconds>(end - start).count() << endl;
   A.output();
   return 0;
}
