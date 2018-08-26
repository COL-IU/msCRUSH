#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <thread>

#include "../class/params.h"
#include "../class/spectrum.h"
#include "../class/core.h"
#include "../class/distance.h"
#include "../utility/io.h"

using namespace Core;
using namespace Utility;
using namespace std;

mutex coutmutex;

struct HyperParams
{
  bool log_normalized;
  bool remove_precursor;
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
  int threads_to_use;
  int window_mz;
  string file_name;
  string result_path;
  string result_prefix;
  string cs_path;
};

void PrintHyperParams(const HyperParams& params) {
  cout << "************ msCRUSH Cluster Params Setting ****************" << endl;
  cout << "cluster iteration: " << params.cluster_iteration << endl;
  cout << "cluster result prefix: " << params.result_prefix << endl;
  cout << "cluster result file path: " << params.result_path << endl;
  cout << "hash dimension: " << params.hash_dimension << endl;
  cout << "hash function num: " << params.hash_func_num << endl;
  cout << "input file name: " << params.file_name << endl;
  cout << "max mz: " << params.max_mz << endl;
  cout << "min mz: " << params.min_mz << endl;
  cout << "min similarity: " << params.min_similarity << endl;
  cout << "mz scale : " << params.mz_scale<< endl;
  cout << "precision: " << params.precision << endl;
  cout << "precursor mass tolerance: " <<  params.mass_tolerance << endl;
  cout << "shared peak mz epsilon: " << params.shared_peak_mz_epsilon << endl;
  cout << "select topk: " << params.select_topk << endl;
  cout << "threds to use: " << params.threads_to_use << endl;
  cout << "window size: " << params.window_mz << endl;
  cout << "********************************************************" << endl;
}

void SetHyperParams(HyperParams* params) {
  (*params).cluster_iteration = 100;  //TODO: Tune this param.
  (*params).hash_func_num = 10;
  (*params).log_normalized = true;
  (*params).mass_tolerance = 0.05;
  (*params).max_mz = 2000;
  (*params).min_mz = 200;
  (*params).min_similarity = 0.65;  //TODO: Tune this param.
  (*params).mz_scale = 1000;
  (*params).precision = 0.8;  //TODO: Tune this param.
  (*params).remove_precursor = false;
  (*params).rounds_num = 4;
  (*params).shared_peak_mz_epsilon = 0.2;
  (*params).select_topk = 5;
  (*params).threads_to_use = 20;
  (*params).window_mz = 100;

  (*params).hash_dimension = int(((*params).max_mz - (*params).min_mz) /
    (*params).precision) + 1;
}

void SaveClusters(const vector<Spectrum* >& spectra_all,
    bool save_spectrum_of_no_charge, string out_file) {
    ofstream out(out_file);
    out << "ID\tTitles" << endl;

    for (int i = 0; i < spectra_all.size(); ++i) {
      if (0 == spectra_all[i]->_charge && !save_spectrum_of_no_charge) {
        continue;
      }
      // Save w/o charge info.
      out << i << "\t" << spectra_all[i]->_component_titles << endl;
    }
    out.close();
}


void p_lsh(vector<Spectrum*> *hash_values, vector<int>* hash_keys,
    hashTable& hash_table, vector<Spectrum*>* unknown_spectra, 
    int hash_func_num, int start_index, int end_index) {
  auto start_time = chrono::high_resolution_clock::now();

  const int buckets = int(pow(2, hash_func_num));
  vector<Spectrum*> local_hash_values;
  vector<int> local_hash_keys;
  
  // Reserve space to speed.
  const int num_points = end_index - start_index;
  local_hash_values.reserve(num_points);
  local_hash_keys.reserve(num_points);

  for (int i = start_index; i < end_index; ++i) {
    auto& spectrum = *(unknown_spectra->at(i));
    int key = LSH::random_projection(spectrum._embeded_peaks, hash_table, false);

    local_hash_values.push_back(unknown_spectra->at(i));
    local_hash_keys.push_back(key);
  }
  swap(*hash_values, local_hash_values);
  swap(*hash_keys, local_hash_keys);
}

void merge_hashtable(vector<vector<Spectrum*>>* final_table, int hash_func_num,
    const vector<vector<Spectrum*>* >& part_hash_values,
    const vector<vector<int>* >& part_hash_keys) {
  const int buckets = int(pow(2, hash_func_num));
  vector<vector<Spectrum*>> local_table(buckets, vector<Spectrum*>());

  for (int i = 0; i < part_hash_keys.size(); ++i) {
    for (int j = 0; j < part_hash_keys[i]->size(); ++j) {
      const auto& key = part_hash_keys[i]->at(j);
      const auto& ptr = part_hash_values[i]->at(j);
      local_table[key].push_back(ptr);
    }
  }
  swap(*final_table, local_table);
}


void merge_spectra(vector<Spectrum*>* final_spectra,
    const vector<vector<Spectrum*>* >& part_spectra) {
  vector<Spectrum*> local_spectra;
  for (const auto& spectra : part_spectra) {
    local_spectra.insert(local_spectra.end(), spectra->begin(), spectra->end());
  }
  swap(*final_spectra, local_spectra);
}

void p_cluster(vector<Spectrum*>* part_spectra, HyperParams* params,
    vector<vector<Spectrum*>>* lsh_table,
    int start_pos, int end_pos, float threshold) {
  
  vector<Spectrum*> local_spectra;
  float distance;
  string component_titles;

  for (int s = start_pos; s < end_pos; ++s) {
    auto& candidates = lsh_table->at(s);

    int size = candidates.size();
    for (int i = 1, j = 0; i < size;) {

      auto& current = *(candidates[i]);
      for (j = 0; j < i; ++j) {

        auto& candidate = *(candidates[j]);

        // Filter using mass_tol and shared_top_peak.
        if (fabs(candidate._precursor_mz - current._precursor_mz) >
          params->mass_tolerance || !current.shareTopPeaks(
            candidate, params->shared_peak_mz_epsilon)) {
          continue;
        }

        distance = Distance::cosine(current._filtered_peaks,
            candidate._filtered_peaks, params->precision);
        if (1 - distance >= threshold) {
          Spectrum* s_new = new Spectrum();
          component_titles =
            current._component_titles + ";" + candidate._component_titles;
          IO::SetConsensus(s_new, current, candidate, params->precision,
              params->select_topk, params->window_mz, params->min_mz,
              params->mz_scale, component_titles, component_titles);
          delete candidates[i];
          delete candidates[j];
          candidates[j] = s_new;
          candidates[i] = candidates[--size];
          candidates[size] = NULL;
          break;
        }
      }
      if (j == i) {
        ++i;
      }
    }
    local_spectra.insert(local_spectra.end(),
        candidates.begin(), candidates.begin() + size);
  }
  swap(*part_spectra, local_spectra);
}

void SplitCommands(int argc, char*argv[], HyperParams* params,
    vector<string>* files) {
  for (int i = 1; i < argc; ++i) {
    if (1 == i) {
      (*params).threads_to_use = stoi(argv[i]);
    } else if (2 == i) {
      (*params).hash_func_num = stoi(argv[i]);
    } else if (3 == i) {
      (*params).cluster_iteration = stoi(argv[i]);
    } else if (4 == i) {
      (*params).min_similarity = stof(argv[i]);
    } else if (5 == i) {
      (*params).result_prefix = string(argv[i]);
    } else {
      (*files).push_back(string(argv[i]));
    }
  }
}

void Cluster(vector<Spectrum*>* unknown_spectra, HyperParams& params, 
    vector<Spectrum*>* spectra_of_no_charge, int charge) {
  auto start_time_total = chrono::high_resolution_clock::now(), 
  start_time = chrono::high_resolution_clock::now(), 
  end_time = chrono::high_resolution_clock::now();

  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>
    (end_time - start_time).count();

  float max_similarity = 0.9;  // Heuristics for maximum cosine similarity.
  float sim_step = (max_similarity - params.min_similarity) / 
    params.cluster_iteration;
  float threshold = max_similarity;
  int cluster_iteration = 0;
  vector<vector<Spectrum*>> lsh_table
    (int(pow(2, params.hash_func_num)), vector<Spectrum*>());

  while (cluster_iteration++ < params.cluster_iteration) {
    threshold -= sim_step;
    cout << "Iteration:\t" << cluster_iteration << ", cos sim threshold:\t" 
      << threshold << endl;

    // Generate LSH.
    hashTable hash_table =
      LSH::generateHashTable(params.hash_func_num, params.hash_dimension);

    // Apply LSH.
    auto start_time_iteration = chrono::high_resolution_clock::now();
    auto start_time_hashing = chrono::high_resolution_clock::now();

    vector<vector<Spectrum*>* > part_hash_values;
    vector<vector<int>* > part_hash_keys;
    for (int i = 0; i < params.threads_to_use; ++i) {
      part_hash_values.push_back(new vector<Spectrum*>());
      part_hash_keys.push_back(new vector<int>());
    }

    int spectra_per_thread = (*unknown_spectra).size() / params.threads_to_use;

    #pragma omp parallel for
    for (int i = 0; i < params.threads_to_use; ++i) {
      int spectra_to_do = spectra_per_thread;
      // If this is the last thread, then take the rest.
      if (i == params.threads_to_use - 1) {
        spectra_to_do = (*unknown_spectra).size() - i * spectra_per_thread;
      }
      int start_index = i * spectra_per_thread;
      p_lsh(part_hash_values[i], part_hash_keys[i], hash_table,
            unknown_spectra, params.hash_func_num, start_index,
            start_index + spectra_to_do);
    }

    // Merge thread results.
    merge_hashtable(&lsh_table, params.hash_func_num, part_hash_values, 
        part_hash_keys);

    for (int i = 0; i < params.threads_to_use; ++i) {
      delete part_hash_values[i];
      delete part_hash_keys[i];
    }
    //part_table.clear();

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_hashing).count();
    cout << "hashing takes secs:\t" << elapsed_read << ", ";

    auto start_time_clustering = chrono::high_resolution_clock::now();
    // Cluster within LSH buckets.
    vector<vector<Spectrum*>*> part_spectra;
    //part_spectra.clear();
    for (int i = 0; i < params.threads_to_use; ++i) {
      part_spectra.push_back(new vector<Spectrum*>());
    }

    int buckets_per_thread = lsh_table.size() / params.threads_to_use;

    #pragma omp parallel for
    for (int i = 0; i < params.threads_to_use; ++i) {

      int start_pos = buckets_per_thread * i;
      int buckets_to_do = buckets_per_thread;

      // If this is the last thread, then take the rest.
      if (i == params.threads_to_use - 1) {
        buckets_to_do = lsh_table.size() - buckets_per_thread * i;
      }

      p_cluster(part_spectra[i], &params, &lsh_table, start_pos,
        start_pos + buckets_to_do, threshold);
    }

    merge_spectra(unknown_spectra, part_spectra);

    for (int i = 0; i < params.threads_to_use; ++i) {
      delete part_spectra[i];
    }

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time_clustering).count();
    cout << "clustering takes secs:\t" << elapsed_read << endl;
    cout << "#spectra after clustering: " << (*unknown_spectra).size() << endl;
    cout << endl;
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time_total).count();
  cout << "msCRUSH algorithm hash+cluster takes (secs): " << elapsed_read << endl;

  // Save clusters.
  start_time = chrono::high_resolution_clock::now();
  cout << "Saving cluster results starts: " << endl;
  //TODO(LEI): customize cluster file name.
  SaveClusters(*unknown_spectra, false, params.result_prefix + "-c" + to_string(charge) + ".txt");
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Save cluster results takes secs: " << elapsed_read << endl;
  
  // Release allocated memory for spectra w/ charge; save pointers to spectra
  // w/o charge to variable 'spectra_of_no_charge'.
  start_time = chrono::high_resolution_clock::now();
  cout << "Releasing memory starts." << endl;
  (*spectra_of_no_charge).clear();
  for (int i = 0; i < (*unknown_spectra).size(); ++i) {
    if (0 == (*(*unknown_spectra)[i])._charge) {
      (*spectra_of_no_charge).push_back((*unknown_spectra)[i]);
    } else {
      delete (*unknown_spectra)[i];
    }
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Releasing memory takes: " << elapsed_read << endl;

  // Report time for all the procedures done above.
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time_total).count();
  cout << "msCRUSH algorithm in total takes (secs): " << elapsed_read << endl;
}

int main (int argc, char *argv[]) {
  if (argc < 7) {
    cout << "Missing parameters, at least 7 params." << endl;
    cout << "Usage: ./mscrush_on_general_charge threads_to_use hash_func_num iteration min_similarity result_prefix mgf_file(s)." << endl;
    return -1;
  }

  HyperParams params;
  SetHyperParams(&params); 

  vector<string> files;
  SplitCommands(argc, argv, &params, &files);

  PrintHyperParams(params);

  auto start_time_total = chrono::high_resolution_clock::now();

  auto start_time = chrono::high_resolution_clock::now();

  bool log_normalized = true, remove_precursor = true;

  vector<Spectrum*> unknown_spectra;
  unordered_map<int, vector<int>> map_spectra_by_charge;
  unordered_map<string, int> map_ms_titles_to_index;

  int unknown_spectra_size = 0;
  cout << "Loading spectra..." << endl;
  for (int i = 0; i < files.size(); ++i) {
    cout << "Loading spectra from file: [" << i+1 <<  "/" << files.size() << "]" << endl;
    cout << "file name: " << files[i] << endl;
    IO::ReadSpectraFromMGF(&unknown_spectra, &map_spectra_by_charge,
        &map_ms_titles_to_index,
        &unknown_spectra_size, files[i], params.mz_scale,
        params.min_mz, params.max_mz, params.precision, params.select_topk,
        params.window_mz, log_normalized, remove_precursor);
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Loading spectra takes secs:\t" << elapsed_read << endl;

  cout << "Now we have #unknown spectra: " << unknown_spectra.size() << endl;

  for (const auto& entry : map_spectra_by_charge) {
    cout << entry.first << " charge: " << entry.second.size() << endl; 
  }

  // Store pointers of spectra w/o charge info.
  vector<Spectrum*> spectra_of_no_charge;
  for (const auto& idx : map_spectra_by_charge[0]) {
    spectra_of_no_charge.push_back(unknown_spectra[idx]);
  }

  // Run clustering on each charge subsequently.
  vector<int> charges{2, 3, 4, 5, 1, 6, 7, 8, 9, 10};
  for (const auto& charge : charges) {
    if (0 == map_spectra_by_charge.count(charge)) {
      cout << "Skip spectra of charge +" << charge << endl;
      continue;
    }
    cout << "Cluster on spectra of charge +" << charge << endl;

    cout << "Load spectra with charge +" << charge << ": ";
    vector<Spectrum*> unknown_spectra_of_interest;
    for (const int& idx : map_spectra_by_charge[charge]) {
      unknown_spectra_of_interest.push_back(unknown_spectra[idx]);
    }
    cout << unknown_spectra_of_interest.size() << endl;

    // Load non-clusterd spectra with no charge info.
    cout << "Load spectra with no charge(+0): ";
    unknown_spectra_of_interest.insert(unknown_spectra_of_interest.end(),
        spectra_of_no_charge.begin(), spectra_of_no_charge.end());
    cout << spectra_of_no_charge.size() << endl;

    cout << "Starting to cluster." << endl;
    Cluster(&unknown_spectra_of_interest, params, &spectra_of_no_charge, 
        charge);
  }

  // Save clusters of spectra w/o charge info.
  cout << "Starting to save spectra w/o charge info." << endl;
  SaveClusters(spectra_of_no_charge, true, params.result_prefix + "-c0.txt");

  // Releasing memory of spectra w/o charge info.
  cout << "Releasing memory of spectra w/o charge info." << endl;
  for (int i = 0; i < spectra_of_no_charge.size(); ++i) {
    delete spectra_of_no_charge[i];
  }

  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time_total).count();
  cout << "In all, msCRUSH takes: " << elapsed_read << endl;
  cout << endl;

  return 0;
}

