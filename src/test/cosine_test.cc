#include <cassert>
#include <chrono>
#include <iostream>
#include <map>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <string>

#include "../class/params.h"
#include "../class/spectrum.h"
#include "../class/core.h"
#include "../class/distance.h"
#include "../utility/io.h"

// #define DEBUG
using namespace Core;
using namespace Utility;
using namespace std;

struct HyperParams
{
  float mass_tolerance;
  float max_mz;
  float min_mz;
  float min_similarity;
  float mz_scale;
  float precision;
  float shared_peak_mz_epsilon;
  int cluster_iteration;
  int hash_func_num;
  int hash_dimension;
  int rounds_num;
  int select_topk;
  int window_mz;
  string file_name;
  string result_path;
};

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 75;  //TODO: Tune this param.
  (*params).mass_tolerance = 0.05;
  (*params).max_mz = 2000;
  (*params).min_mz = 200;
  (*params).min_similarity = 0.55;  //TODO: Tune this param.
  (*params).mz_scale = 1000;
  (*params).precision = 0.8;  //TODO: Tune this param.
  (*params).rounds_num = 4;
  (*params).shared_peak_mz_epsilon = 0.2;
  (*params).hash_func_num = 12;
  (*params).hash_dimension = int(((*params).max_mz - (*params).min_mz) /
    (*params).precision) + 1;
  (*params).select_topk = 5;
  (*params).window_mz = 100;
}

int main (int argc, char *argv[]) {
  HyperParams params;
  SetHyperParams(&params); 
  params.file_name = argv[1];
  cout << "mgf file name: " << params.file_name << endl;

  auto start_time = chrono::high_resolution_clock::now();
  cout << "Loading MS2 data from disk." << endl;
  const vector<Spectrum> &unknown_spectra_all =
    IO::ReadSpectraFromMGF(params.file_name, params.mz_scale, params.min_mz,
        params.max_mz, params.precision, params.select_topk, params.window_mz);
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Loading spectra takes secs:\t" << elapsed_read << endl;


  double max_diff = 0;
  int ii = -1, jj = -1;
  for (int i = 0; i < unknown_spectra_all.size(); ++i) {
    for (int j = i + 1; j < unknown_spectra_all.size(); ++j) {
      auto s1 = unknown_spectra_all[i];
      // cout << s1 << endl;
      auto s2 = unknown_spectra_all[j];
      if (fabs(s1._precursor_mz - s2._precursor_mz) > params.mass_tolerance) {
        continue;
      }

      // cout << s2 << endl;
      cout << s1._title << endl << s2._title << endl;
      cout << "old cosine dist: " << endl;
      double old_dist = Distance::cosine(s1._embeded_peaks, s2._embeded_peaks);
      cout << old_dist << endl;
      cout << "new cosine dist: "<< endl;
      double new_dist = Distance::cosine(s1._filtered_peaks, s2._filtered_peaks, params.precision);
      cout << new_dist << endl;
      //assert(new_dist <= old_dist);
      double temp = fabs(old_dist - new_dist);
      if (temp > max_diff) {
        max_diff = temp;
        ii = i;
        jj = j;
      }
      cout << "diff: " << temp << endl;
      cout <<  "-----------" << endl;
    }
  }
  cout << "max diff: " << max_diff << endl;
  cout << unknown_spectra_all[ii]._title << endl;
  cout << unknown_spectra_all[jj]._title<< endl;

  return 0;
}
  
