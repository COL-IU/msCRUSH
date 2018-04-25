#include <iostream>
#include "../class/peak.h"
using namespace std;
using namespace Core;

void print(const Peak& p) {
  cout << p._mz << endl;
  cout << p._intensity << endl;
  cout << p._count << endl;
}

void print(const EmbededPeak& p) {
  cout << p._idx << endl;
  cout << p._intensity << endl;
  cout << p._count << endl;
}

int main () {
  EmbededPeak p;
  print(p);

  EmbededPeak pp(1, 10, 1);
  print(pp);
  pp = p;
  print(pp);
  return 0;
}
