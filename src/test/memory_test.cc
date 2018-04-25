#include <iostream>
#include "../class/spectrum.h"

using namespace std;
using namespace Core;

vector<vector<Spectrum>* > indexed_spectra;

int main (int argc, char *argv[]) {
  int num_limit = 100000;

  for (int i = 0; i < 100; ++i) {
    indexed_spectra.push_back(new vector<Spectrum>(num_limit));
    cout << "allocating memory: #" << i + 1 << endl;
  }

  return 0;
}
