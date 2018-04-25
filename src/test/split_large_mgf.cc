#include <iostream>

#include "../utility/io.h"
#include "../class/distance.h"
using namespace std;
using namespace Core;
using namespace Utility;

void WriteCS(const vector<Spectrum* >& spectra_all, string cs_path) {
    string base_path = cs_path;
    ofstream out;
    int max_spectra_per_file = 400000;
    for (int i = 0; i < spectra_all.size(); ++i) {

      if (i % max_spectra_per_file == 0) {
        out.close();
        out.clear();
        out.open(base_path + to_string(i / max_spectra_per_file) + ".mgf");
      }
      const auto& spectrum = *(spectra_all[i]);
      out << "BEGIN IONS" << endl;
      out << "PEPMASS=" << spectrum._precursor_mz << endl;
      out << "CHARGE=" << spectrum._charge <<"+" << endl;
      out << "TITLE=" << spectrum._title << endl;

      for (const auto& peak : spectrum._filtered_peaks) {
        out << peak._mz << "\t" << peak._intensity << endl;
      }
      out << "END IONS" << endl;
      out << endl;
    }
    out.close();
}

int main (int argc, char *argv[]) {
  string file_name = argv[1];
  string cs_path = argv[2];
  cout << "file name: " << file_name << endl;

  vector<Spectrum*> indexed_spectra;
  float min_mz = 200, max_mz = 2000, precision = 0.8, scale = 1000;
  int select_topk = 5, window_mz = 100;
  int size = 0;
  IO::ReadSpectraFromMGF(&indexed_spectra, &size, file_name, scale, min_mz,
      max_mz, precision, select_topk, window_mz);


  cout << "size: " << size << endl;
  cout << "#read spectra: " << indexed_spectra.size() << endl;

  WriteCS(indexed_spectra, cs_path);

  //cout << *(indexed_spectra.front()) << endl;
  //cout << *(indexed_spectra.back()) << endl;


  //Spectrum tmp;
  //IO::SetConsensus(&tmp, indexed_spectra[0], indexed_spectra[1],
  //    precision, select_topk, window_mz, min_mz, scale, "tmp", "tmp");
  //indexed_spectra.push_back(tmp);

  //Spectrum tmp;
  //IO::SetConsensus(&tmp, indexed_spectra[0], indexed_spectra[2],
  //    precision, select_topk, window_mz, min_mz, scale, "tmp", "tmp");

  //Spectrum tmp2;
  //IO::SetConsensus(&tmp2, tmp, indexed_spectra[1],
  //    precision, select_topk, window_mz, min_mz, scale, "consensus", "consensus");

  //indexed_spectra.push_back(tmp);
  //indexed_spectra.push_back(tmp2);

  //for (const auto& spectrum : indexed_spectra) {
  //  cout << "BEGIN IONS" << endl;
  //  cout << "PEPMASS=" << spectrum._precursor_mz << endl;
  //  cout << "CHARGE=" << spectrum._charge <<"+" << endl;
  //  cout << "TITLE=" << spectrum._title << endl;
  //  
  //  for (const auto& peak : spectrum._filtered_peaks) {
  //    cout << peak._mz << "\t" << peak._intensity << endl;
  //  }
  //  cout << "END IONS" << endl;
  //  cout << endl;
  //}

  //cout << "cosine similarity." << endl;
  //for (int i = 0; i < indexed_spectra.size(); ++i) {
  //  for (int j = i + 1; j < indexed_spectra.size(); ++j) {
  //    auto& si = *(indexed_spectra[i]);
  //    auto& sj = *(indexed_spectra[j]);
  //    cout << si._title << " vs " << sj._title << ": ";
  //    cout << 1 - Distance::cosine(si._filtered_peaks, sj._filtered_peaks, precision) << endl;
  //  }
  //}
  
  
  for (int i = 0; i < size; ++i) {
    delete indexed_spectra[i];
  }
  return 0;
}
