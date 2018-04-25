#include <iostream>
#include <vector>
using namespace std;

void worker (vector<int*> &v) {
  for (int i = 0; i < 5; ++i) {
    int* s = new int(i);
    v.push_back(s);
  }
}
int main () {
  vector<int*> v;
  worker(v);
  for (int i = 0; i < 5; ++i) {
    cout << v[i] << " : " <<  *(v[i]) << endl;
    delete v[i];
  }

  cout << v.size() << endl;
  return 0;
}

