#include <iostream>

#include "../utility/io.h"
using namespace std;
using namespace Core;
using namespace Utility;

int main (int argc, char *argv[]) {
  string file_name = argv[1];
  //string file_name = "../../data/Adult_Liver_8times.mgf";
  //string file_name = "../../data/Adult_Liver.mgf";
  cout << "file name: " << file_name << endl;
  float min_mz = 0, max_mz = 99999, precision = 0.8, scale = 1000;
  int select_topk = 5, window_mz = 100;
  auto indexed_spectra = IO::ReadSpectraFromMGF(file_name, scale, min_mz, max_mz, precision, select_topk, window_mz);
  cout << "#read spectra: " << indexed_spectra.size() << endl;

  for (const auto& spectrum : indexed_spectra) {
    cout << spectrum << endl;
  }
  Spectrum consensus;
  IO::SetConsensus(&consensus, indexed_spectra, precision, select_topk, window_mz, min_mz, scale);

  for (const auto& peak : consensus._embeded_peaks) {
    cout << peak._idx * precision << "\t" << peak._intensity << endl;
  }
  return 0;
}
