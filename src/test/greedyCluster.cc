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
  (*params).cluster_iteration = 50;  //TODO: Tune this param.
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

  (*params).file_name = "/home/wang558/2_Adult_Liver.mgf";
  (*params).result_path = "/home/wang558/toy";

}


unordered_map<int, Spectrum> unknown_spectra;

int main (int argc, char *argv[]) {
  HyperParams params;
  SetHyperParams(&params); 

  //params.file_name = argv[1];

  auto start_time = chrono::high_resolution_clock::now();
  cout << "Loading MS2 data from disk." << endl;
  //vector<Spectrum> unknown_spectra;
  IO::ReadSpectraFromMGF(&unknown_spectra, params.file_name, params.mz_scale,
      params.min_mz, params.max_mz, params.precision, params.select_topk,
      params.window_mz);

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Loading spectra takes secs:\t" << elapsed_read << endl;

  return 0;
}

