#include <chrono>
#include <iostream>

#include "../utility/io.h"
#include "../class/distance.h"
using namespace std;
using namespace Core;
using namespace Utility;

//void WriteConsensusSpectrum(const Spectrum& spectrum, ofstream* ofs) {
void WriteConsensusSpectrum(const Spectrum& spectrum, string title_prefix, 
    int index, int charge, ofstream* ofs) {
  (*ofs) << "BEGIN IONS" << endl;
  (*ofs) << "PEPMASS=" << spectrum._precursor_mz << endl;
  if (0 != spectrum._charge && -1 != spectrum._charge) {
    (*ofs) << "CHARGE=" << spectrum._charge <<"+" << endl;
  }
  //(*ofs) << "TITLE=" << spectrum._title << endl;
  (*ofs) << "TITLE=" << title_prefix << "." << to_string(index) << "." <<
    to_string(index) << "." << to_string(charge) << ".dta" << endl;
  for (const auto& peak : spectrum._filtered_peaks) {
    (*ofs) << peak._mz << "\t" << peak._intensity << endl;
  }
  (*ofs) << "END IONS" << endl;
  (*ofs) << endl;
}

void WriteConsensusSpectra(const vector<Spectrum>& spectra, string title_prefix,
    int charge, ofstream* ofs) {
    for (int i = 0; i < spectra.size(); ++i) {
      if (i && i % 10000 == 0) {
        cout << i << " write done." << endl;
      }
      const auto& spectrum = (spectra[i]);
      //WriteConsensusSpectrum(spectrum, ofs);
      WriteConsensusSpectrum(spectrum, title_prefix, i, charge, ofs);
    }
}
void SplitCommands(int argc, char*argv[], string* title_prefix,
    string* cs_path_prefix, vector<string>* cluster_files, 
    vector<string>* mgf_files) {
  for (int i = 1; i < argc; ++i) {
    if (1 == i) {
      *title_prefix = string(argv[i]);
    } else if (2 == i) {
      *cs_path_prefix = string(argv[i]);
    } else if (string(argv[i]).rfind(".txt") != string::npos || 
        string(argv[i]).rfind(".sc") != string::npos){
      (*cluster_files).push_back(string(argv[i]));
    } else if (string(argv[i]).rfind(".mgf") != string::npos){
      (*mgf_files).push_back(string(argv[i]));
    } else {
      cout << "Errors splitting commands!" << endl;
      exit(-1);
    }
  }
}

int main (int argc, char *argv[]) {
  if (argc < 5) {
    cout << "Missing parameters, at least 5 params." << endl;
    cout << "Usage: ./program_exe cs_title_prefix cs_path_prefix splash_clusters_name(s) file(s)." << endl;
    return -1;
  }
  auto start_time_total = chrono::high_resolution_clock::now();

  auto start_time = chrono::high_resolution_clock::now();

  int charge;
  string cs_path_prefix, title_prefix;
  vector<string> mgf_files, cluster_files;
  SplitCommands(argc, argv, &title_prefix, &cs_path_prefix, &cluster_files,
      &mgf_files);

  vector<Spectrum*> indexed_spectra;
  unordered_map<int, vector<int>> map_spectra_by_charge;
  unordered_map<string, int> map_ms_titles_to_index;

  float min_mz = 0, max_mz = 2000, precision = 0.05, scale = 1000;
  bool log_normalized = false, remove_precursor = true;
  int select_topk = 30, window_mz = 100;
  int size = 0;

  for (int i = 0; i < mgf_files.size(); ++i) {
    cout << "Loading spectra from mgf file: [" << i+1 <<  "/" << mgf_files.size() << "]" << endl;
    cout << "file name: " << mgf_files[i] << endl;
    IO::ReadSpectraFromMGF(&indexed_spectra, &map_spectra_by_charge,
        &map_ms_titles_to_index, &size, mgf_files[i], scale,
        min_mz, max_mz, precision, select_topk,
        window_mz, log_normalized, remove_precursor);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Loading spectra takes sec:\t" << elapsed_read << endl;

  //cout << "size: " << size << endl;
  cout << "#read spectra: " << indexed_spectra.size() << endl;

  // Read clusters from file.
  start_time = chrono::high_resolution_clock::now();
  int file_cnt = 0;
  for (const string cluster_file : cluster_files) {
    cout << "Generating CS for cluster file: [" << ++file_cnt <<  "/" << cluster_files.size() << "]" << endl;
    cout << "file name:\t" << cluster_file << endl;
    string tmp = cluster_file.substr(0, cluster_file.rfind("."));
    charge = stoi(tmp.substr(tmp.find_last_of('c') + 1));

    ifstream inf(cluster_file);
    ofstream ofs(cs_path_prefix + to_string(charge) + ".cs.mgf");
    string line, delimiter = ";";
    int cnt = 0;

    // Skip header.
    getline(inf, line);

    while (getline(inf, line)) {
      std::string token;
      istringstream iss(line.substr(line.find_first_of('\t') + 1));

      Spectrum s_new;
      if (cnt + 1 % 100000 == 0) {
        cout << cnt + 1 << " clustering done." << endl;
      }
      while (getline(iss, token, ';')) {
        auto& candidate = *(indexed_spectra[map_ms_titles_to_index[token]]);
        if (s_new._title == "default") {
          s_new = candidate;
        } else {
          string component_titles =
            s_new._component_titles + ";" + candidate._component_titles;
          IO::SetConsensus(&s_new, candidate, s_new, precision, select_topk,
              window_mz, min_mz, scale, component_titles, component_titles);
        }
      }
      WriteConsensusSpectrum(s_new, title_prefix + to_string(charge), 
          cnt, charge, &ofs);
      ++cnt;
    }
    cout << cnt << " clustering done." << endl;

    ofs.close();
    inf.close();
  }

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();

  cout << "Writing consensus spectra done with secs:\t" << elapsed_read << endl;

  cout << "Starting to release memory." << endl;
  
  start_time = chrono::high_resolution_clock::now();
  for (int i = 0; i < size; ++i) {
    delete indexed_spectra[i];
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Completed releasing memory with secs:\t" << elapsed_read << endl;
  return 0;
}
