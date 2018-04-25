#include <iostream>
#include <stdio.h>
using namespace std;

int main () {
  float a, b;
  //string line = "1.3  4.5";
  string line = "end ions";
  int n = sscanf(line.c_str(), "%f%f", &a, &b);
  cout << "return: " << n << endl;
  cout << a << endl;
  cout << b << endl;

  return 0;
}
