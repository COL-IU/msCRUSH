#include <iostream>

#include "../utility/io.h"
#include "../class/distance.h"
using namespace std;
using namespace Core;
using namespace Utility;

void WriteCS(const vector<Spectrum* >& spectra_all, string cs_path) {
    string out_file = cs_path;
    ofstream out(out_file);
    for (int i = 0; i < spectra_all.size(); ++i) {
      const auto& spectrum = *(spectra_all[i]);
      out << "BEGIN IONS" << endl;
      out << "PEPMASS=" << spectrum._precursor_mz << endl;
      out << "CHARGE=" << spectrum._charge <<"+" << endl;
      out << "TITLE=" << spectrum._title << endl;

      for (const auto& peak : spectrum._filtered_peaks) {
        //TODO: intention is to reverse naturl logarithm done in preprocessing.
        out << peak._mz << "\t" << peak._intensity << endl;
        //out << peak._mz << "\t" << exp(peak._intensity/100.) << endl;
      }
      out << "END IONS" << endl;
      out << endl;
    }
    out.close();
}

int main (int argc, char *argv[]) {
  string file_name = argv[1];
  string normalized = string(argv[2]);

  if (normalized != "yes" && normalized != "no") {
    cout << "you have to input 'yes' or 'no' to indicate log_normalized!" << endl;
    exit(0);
  }
  bool log_normalized = normalized == "yes";

  vector<Spectrum*> indexed_spectra;
  float min_mz = 0, max_mz = 2000, precision = 0.05, scale = 1000;
  int select_topk = 30, window_mz = 100;
  int size = 0;
  IO::ReadSpectraFromMGF(&indexed_spectra, &size, file_name, scale, min_mz,
      max_mz, precision, select_topk, window_mz, log_normalized, false);

  cout << "file name: " << file_name << endl;
  cout << "log normalized ? " << (log_normalized ? "True": "False" ) << endl;
  cout << "size: " << size << endl;
  cout << "#read spectra: " << indexed_spectra.size() << endl;

  // cout << *(indexed_spectra.back()) << endl;
  // cout << indexed_spectra.back()->_title.length() << endl;

  //WriteCS(indexed_spectra, cs_path);

  //cout << *(indexed_spectra.front()) << endl;
  //cout << *(indexed_spectra.back()) << endl;


  //Spectrum* tmp  = new Spectrum();
  //IO::SetConsensus(tmp, *(indexed_spectra[0]), *(indexed_spectra[1]),
  //    precision, select_topk, window_mz, min_mz, scale, "consensus", "consensus");
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
  //  cout << "PEPMASS=" << (*spectrum)._precursor_mz << endl;
  //  cout << "CHARGE=" << (*spectrum)._charge <<"+" << endl;
  //  cout << "TITLE=" << (*spectrum)._title << endl;
  //  
  //  for (const auto& peak : (*spectrum)._filtered_peaks) {
  //    cout << peak._mz << "\t" << peak._intensity << endl;
  //  }
  //  cout << "END IONS" << endl;
  //  cout << endl;
  //}

  cout << "cosine similarity." << endl;
  unordered_map<string, string> mymap = {
    {"Adult_Monocytes_bRP_Velos_31_f19.1339.1339.2.dta", "SASPSSVEVHPVLEK"},
    {"Adult_CD4Tcells_bRP_Velos_29_f14.1955.1955.2.dta", "WLHNEDQMAVEK"},
    {"Adult_CD4Tcells_bRP_Velos_29_f15.1728.1728.2.dta", "INYDLPNNR"},
    {"Adult_CD4Tcells_bRP_Velos_29_f01.2044.2044.2.dta", "AAPPPPPPPPPLESSPR"},
    {"Adult_Monocytes_bRP_Velos_31_f17.1462.1462.2.dta", "GGGGNFGPGPGSNFR"},
    {"Adult_Lung_Gel_Velos_13_f08.4135.4135.2.dta", "WLHGLESQIQSDDYGK"},

  };

  for (int i = 0; i < indexed_spectra.size(); ++i) {
    for (int j = i + 1; j < indexed_spectra.size(); ++j) {
      auto& si = *(indexed_spectra[i]);
      auto& sj = *(indexed_spectra[j]);
      //cout << si._title << " vs " << sj._title << ": ";
      cout << mymap[si._title] << " vs " << mymap[sj._title] << ": ";
      cout << 1 - Distance::cosine(si._filtered_peaks, sj._filtered_peaks, precision) << endl;
    }
  }

  //WriteCS(indexed_spectra, "tmp-nolog.mgf");
  
  
  for (int i = 0; i < size; ++i) {
    delete indexed_spectra[i];
  }
  return 0;
}
