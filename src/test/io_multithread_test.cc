#include <iostream>

#include "../utility/io.h"
using namespace std;
using namespace Core;
using namespace Utility;

int main (int argc, char *argv[]) {
  // string file_name = "../../data/Adult_Liver.mgf";
  string file_name = argv[1];
  int spectrum_per_thread = stoi(argv[2]);
  cout << "file name: " << file_name << endl;
  float min_mz = 0, max_mz = 99999, precision = 0.02, scale = 1000;
  int select_topk = 5, window_mz = 100;

  auto indexed_spectra = IO::ReadSpectraFromMGF(-1, spectrum_per_thread, file_name, scale, min_mz, max_mz, precision, select_topk, window_mz);

  cout << "#all spectra: " << indexed_spectra.size() << endl; 
 
  return 0;
}
