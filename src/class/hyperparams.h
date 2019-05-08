#ifndef __CORE_HYPERPARAMS_H__
#define __CORE_HYPERPARAMS_H__

#include <iostream>
#include <iomanip>
#include <string>

using std::cout;
using std::endl;
using std::string;
using std::ostream;

namespace Core {
class HyperParams{
 public:
  bool log_normalized;
  bool remove_precursor;
  float precursor_mass_tolerance;
  float max_mz;
  float min_mz;
  float min_similarity;
  float mz_scale;
  float precision;
  float shared_peak_mz_epsilon;
  int hash_func_num;
  int hash_dimension;
  int iteration;
  int select_topk;
  int threads_to_use;
  int window_mz;
  string consensus_file_prefix;
  string clustering_file_prefix;

  HyperParams() {
    consensus_file_prefix = "consensus";
    clustering_file_prefix = "clustering";
    hash_func_num = 10;
    iteration = 100;
    log_normalized = true;
    max_mz = 2000;
    min_mz = 200;                                                                                                  
    min_similarity = 0.65;
    mz_scale = 1000;
    precision = 0.8;
    precursor_mass_tolerance = 0.05;
    remove_precursor = true;
    shared_peak_mz_epsilon = 0.2;
    select_topk = 5;
    threads_to_use = 20;
    window_mz = 100;
    hash_dimension = int((max_mz - min_mz) / precision) + 1;
  }

  friend ostream& operator<<(ostream& os, const HyperParams& params);
};

inline ostream& operator<<(ostream& os, const HyperParams& params) {
  os << "--Output HyperParams Begins--" << endl;
  os << "<bool> peak intensity log normalized: "
      << (params.log_normalized ? "Yes" : "No") << endl;
  os << "<bool> remove precursor peak: "
      << (params.remove_precursor ? "Yes" : "No") << endl;
  os << "<float> precursor mass tolerance(Da): "
      << params.precursor_mass_tolerance << endl;
  os << "<float> max peak mz: " << params.max_mz << endl;
  os << "<float> min peak mz: " << params.min_mz << endl;
  os << "<float> min cosine sim for clustering: " << params.min_similarity << endl;
  os << "<float> mz scale after nomrlization: " << params.mz_scale << endl;
  os << "<float> fragment precision: " << params.precision << endl;
  os << "<float> shared peak mz eps: " << params.shared_peak_mz_epsilon << endl;
  os << "<int> clustering iteration: " << params.iteration << endl;
  os << "<int> hash functions per hash table: " << params.hash_func_num << endl;
  os << "<int> hash function dimension: " << params.hash_dimension << endl;
  os << "<int> topK strongest peaks per bin: " << params.select_topk << endl;
  os << "<int> threads to use: " << params.threads_to_use << endl;
  os << "<int> window(bin) size: " << params.window_mz << endl;
  os << "<string> consensus spectra file prefix: "
      << params.consensus_file_prefix << endl;
  os << "<string> clustering result file prefix: "
      << params.clustering_file_prefix << endl;
  os << "--Output HyperParams Ends--" << endl;
  return os;
}

}  // namespace Core
#endif
