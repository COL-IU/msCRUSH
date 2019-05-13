#include <chrono>
#include <iostream>

#include "../class/distance.h"
#include "../class/hyperparams.h"
#include "../utility/cmdparser.h"
#include "../utility/io.h"
using namespace std;
using namespace Core;
using namespace Utility;

void WriteConsensusSpectrum(const Spectrum& spectrum, string title_prefix, 
    int index, int charge, ofstream* ofs) {
  (*ofs) << "BEGIN IONS" << "\n";
  (*ofs) << "PEPMASS=" << spectrum._precursor_mz << "\n";
  if (0 != spectrum._charge && -1 != spectrum._charge) {
    (*ofs) << "CHARGE=" << spectrum._charge <<"+" << "\n";
  }

  (*ofs) << "TITLE=" << title_prefix << "." << to_string(index) << "." <<
    to_string(index) << "." << to_string(charge) << ".dta" << "\n";
  for (const auto& peak : spectrum._filtered_peaks) {
    (*ofs) << peak._mz << "\t" << peak._intensity << "\n";
  }
  (*ofs) << "END IONS" << "\n";
  (*ofs) << "\n";
}

void WriteConsensusSpectra(const vector<Spectrum>& spectra, string title_prefix,
    int charge, ofstream* ofs) {
    for (int i = 0; i < spectra.size(); ++i) {
      if (i && i % 10000 == 0) {
        cout << i << " write done." << "\n";
      }
      const auto& spectrum = (spectra[i]);
      WriteConsensusSpectrum(spectrum, title_prefix, i, charge, ofs);
    }
}

void ParseCommands(
    const cli::Parser& parser, 
    int *d, 
    char* separator,
    string* title_prefix,
    string* file_prefix, 
    vector<string>* cluster_files,
    vector<string>* mgf_files) {
  *d = parser.get<int>("d");
  *separator = parser.get<string>("s")[0];
  *title_prefix = parser.get<string>("t");
  *file_prefix = parser.get<string>("p");
  *cluster_files = parser.get<vector<string> >("c");
  *mgf_files = parser.get<vector<string> >("f");
}

void configure_parser(cli::Parser& parser) {  
  parser.set_optional<int>("d", "decimal", 3, "Decimal places for numbers.");
  parser.set_optional<string>("p", "consensus_path_prefix", "consensus",
                              "Consensus result file prefix to write into");
  parser.set_optional<string>("s", "separator", "|", 
                            "Delimiter to separate MS2 titles in clusters");
  parser.set_optional<string>("t", "consensus_title", "CONSENSUS", 
                              "Consensus spectrum title prefix.");

  parser.set_required<vector<string> >("c", "clusters", 
                                      "Clustering files by msCRUSH.");
  parser.set_required<vector<string> >("f", "files", "MGF files.");
}

void ConfigureHyperParams(HyperParams* params) {
  params->min_mz = 0;
  params->max_mz = 2000;
  params->precision = 0.05;
  params->mz_scale = 1000;
  params->log_normalized = false;
  params->remove_precursor = true;
  params->select_topk = 30;
  params->window_mz = 100;
}

int main (int argc, char *argv[]) {
  cli::Parser parser(argc, argv);
  configure_parser(parser);
  parser.run_and_exit_if_error();

  auto start_time_total = chrono::high_resolution_clock::now();
  auto start_time = chrono::high_resolution_clock::now();

  int charge, decimal_place;
  char separator;
  string file_prefix, title_prefix;
  vector<string> mgf_files, cluster_files;
  ParseCommands(parser, &decimal_place, &separator, &title_prefix, &file_prefix, 
                &cluster_files, &mgf_files);

  cout << fixed << showpoint;
  cout << setprecision(decimal_place);
  cout << "Delimiter to separate MS2 titles in clusters: " << separator << endl;

  vector<Spectrum*> indexed_spectra;
  unordered_map<int, vector<int>> map_spectra_by_charge;
  unordered_map<string, int> map_ms_titles_to_index;

  HyperParams params;
  ConfigureHyperParams(&params);

  int size = 0;

  for (int i = 0; i < mgf_files.size(); ++i) {
    cout << "Loading spectra from mgf file: [" << i+1 <<  "/" << mgf_files.size() << "]" << "\n";
    cout << "file name: " << mgf_files[i] << "\n";
    IO::ReadSpectraFromMGF(
        &indexed_spectra, &map_spectra_by_charge, &map_ms_titles_to_index, 
        &size, mgf_files[i], params.mz_scale, params.min_mz, params.max_mz,
        params.precision, params.select_topk, params.window_mz, 
        params.log_normalized, params.remove_precursor);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Loading spectra takes sec:\t" << elapsed_read << "\n";

  //cout << "size: " << size << "\n";
  cout << "#read spectra: " << indexed_spectra.size() << "\n";

  // Read clusters from file.
  start_time = chrono::high_resolution_clock::now();
  int file_cnt = 0;
  for (const string cluster_file : cluster_files) {
    cout << "Generating Consensus Spectra for cluster file: [" 
        << ++file_cnt <<  "/" << cluster_files.size() << "]" << "\n";
    cout << "Clustering file name:\t" << cluster_file << "\n";
    string tmp = cluster_file.substr(0, cluster_file.rfind("."));
    charge = stoi(tmp.substr(tmp.find_last_of('c') + 1));

    ifstream inf(cluster_file);
    ofstream ofs(file_prefix + "-c" + to_string(charge) + ".mgf");
    string line;
    int cnt = 0;

    // Skip header.
    getline(inf, line);

    while (getline(inf, line)) {
      std::string token;
      istringstream iss(line.substr(line.find_first_of('\t') + 1));

      Spectrum s_new;
      if (cnt + 1 % 100000 == 0) {
        cout << cnt + 1 << " clustering done." << "\n";
      }
      while (getline(iss, token, separator)) {
        auto& candidate = *(indexed_spectra[map_ms_titles_to_index[token]]);
        if (s_new._title == "NA") {
          s_new = candidate;
        } else {
          string component_titles =
            s_new._component_titles + separator + candidate._component_titles;
          IO::SetConsensus(
              &s_new, candidate, s_new, params.precision, params.select_topk,
              params.window_mz, params.min_mz, params.mz_scale,
              component_titles, component_titles);
        }
      }
      WriteConsensusSpectrum(s_new, title_prefix, cnt, charge, &ofs);
      ++cnt;
    }
    cout << cnt << " clustering done." << "\n";

    ofs.close();
    inf.close();
  }

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();

  cout << "Writing consensus spectra done with secs:\t" << elapsed_read << "\n";

  cout << "Starting to release memory." << "\n";
  
  start_time = chrono::high_resolution_clock::now();
  for (int i = 0; i < size; ++i) {
    delete indexed_spectra[i];
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Completed releasing memory with secs:\t" << elapsed_read << "\n";
  return 0;
}
