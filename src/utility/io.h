#ifndef __UTILITY_IO_H__
#define __UTILITY_IO_H__

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fcntl.h>
#include <fstream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unistd.h>
#include <utility>
#include <vector>

#include "../class/spectrum.h"
#include "../class/peak.h"
#include "params.h"

using namespace std;
using namespace Core;

namespace Utility {
class IO {
 public:
  // Use customized char* to float converter, it can support up to precision
  // 0.00001. Changes to [strtof] might be needed to handle higher precision.
  static void ReadSpectraFromMGF(vector<Spectrum*>* indexed_spectra, 
      unordered_map<int, vector<int>>* map_spectra_by_charge,
      unordered_map<string, int>* map_ms_titles_to_index,
      int* spectra_size, string file_name, float scale, float min_mz,
      float max_mz, float precision, int select_topk, int window_mz,
      bool peak_normalized, bool remove_precursor, bool verbose=false);

  // [Deprecated].
  static void ProcessSpectra(vector<Spectrum>* spectra, int start, int end,
      float scale, float min_mz, float precision, int select_topk,
      float window_mz);

  // [Deprecated].
  static void ProcessSpectra(const vector<vector<string>>& chunks, float scale,
      float min_mz, float precision, int select_topk, float window_mz);

  static void SetSpectrum(Spectrum* spectrum,
      bool is_clustered, bool is_consensus, int charge,
      int count, const EmbededPeaks& embeded_peaks, const Peaks& raw_peaks,
      const Peaks& filtered_peaks,
      string peptide, float pre_mz, string title, string component_titles,
      const vector<float>& top_peak_mz);

  // Normalize peak intensity with the maximum set to 'scale'.
  static void Normalize(Peaks* peaks, float scale);

  // Remove adjacent peaks within mz tolerance, keep the strongest peak.
  static void RemoveAdjacentPeaks(Peaks* peaks, float mz_tolerance);

  // Embed peaks.
  static void Embed(EmbededPeaks* embeded_peaks, const Peaks& peaks,
      float min_mz, float precision, float scale);
  
  // Set consensus spectrum, using just two spectra.
  static void SetConsensus(Spectrum* consensus, const Spectrum& s1,
      const Spectrum& s2, float precision, int topK, float bin_size,
      float min_mz, float scale, string title, string component_titles);

  static void MergeTwoPeaks(const Peaks& p1, const Peaks& p2, Peaks* peaks);

  // Adjust intensity according to frequency.
  static void AdaptPeakIntensities(Peaks* peaks, int nSpectra);

  static void BinTopKPeak(Peaks* top_peaks, const Peaks& peaks, int peaks_size,
      int topK, float bin_size);

  // [Deprecated].
  //static vector<float> SelectTopPeakMZ(
  //    const EmbededPeaks& peaks, float precision,int topK = 5);

  static vector<float> SelectTopPeakMZ(const Peaks& peaks, int topK = 5);
  
  // Customized char* to float converter, only work with precision up to 0.00001
  static float convert(char const* source, char ** endPtr);
};

}  // namespace Utility
#endif
