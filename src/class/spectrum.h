#ifndef __CORE_SPECTRUM_H__
#define __CORE_SPECTRUM_H__

#include <cmath>
#include <iostream>
#include <string>

#include "params.h"

using std::endl;
using std::fabs;
using std::ostream;
using std::string;

namespace Core {
class Spectrum {
 public:
  bool _is_clustered;
  bool _is_consensus;
  float _precursor_mz;
  int _charge;
  // Number of spectra contributing to this spectrum.
  int _count;
  // Ascendingly sorted mz of topK strongest peaks.
  vector<float> _top_peak_mz;

  // Store the titles of its component spectra.
  // For instance:
  //  1. If this spectrum is not a consensus spectrum,
  //    then _title = _component_titles.
  //  2. If this spectrum is a consensus spectrum composed of s_a, s_b, s_c,
  //    then _component_titles =
  //    ";".join(_title of s_a, _title of s_b, _title of s_c).
  string _component_titles;
  string _peptide;
  string _title;
  // Embedded peaks, with mz-to-idx and intensity.
  EmbededPeaks _embeded_peaks;
  // Filtered peaks, such as Top5 per 100 Da.
  Peaks _filtered_peaks;
  // Raw peaks, with mz in [min_mz, max_mz].
  Peaks _raw_peaks;
  
  Spectrum(): _is_clustered(false), _is_consensus(false), _precursor_mz(0),
  _charge(0), _count(1), _peptide("default"), _title("default"),
  _component_titles("default") {};

  ~Spectrum() {}

  Spectrum(const Spectrum& sp) {
    _charge = sp._charge;
    _component_titles = sp._component_titles;
    _count = sp._count;
    _embeded_peaks = sp._embeded_peaks;
    _filtered_peaks = sp._filtered_peaks;
    _is_clustered = sp._is_clustered;
    _is_consensus = sp._is_consensus;
    _peptide = sp._peptide;
    _precursor_mz = sp._precursor_mz;
    _raw_peaks = sp._raw_peaks;
    _title = sp._title;
    _top_peak_mz = sp._top_peak_mz;
  }

  Spectrum& operator=(const Spectrum& sp) {
    _charge = sp._charge;
    _component_titles = sp._component_titles;
    _count = sp._count;
    _embeded_peaks = sp._embeded_peaks;
    _filtered_peaks = sp._filtered_peaks;
    _is_clustered = sp._is_clustered;
    _is_consensus = sp._is_consensus;
    _peptide = sp._peptide;
    _precursor_mz = sp._precursor_mz;
    _raw_peaks = sp._raw_peaks;
    _title = sp._title;
    _top_peak_mz = sp._top_peak_mz;
    return *this;
  }
  
  friend ostream& operator<<(ostream& os, const Spectrum& spectrum);

  bool shareTopPeaks(const Spectrum& other, float epsilon) {
    int i = 0, j = 0;
    while (i < _top_peak_mz.size() && j < other._top_peak_mz.size()) {
      if (fabs(_top_peak_mz[i] - other._top_peak_mz[j]) < epsilon) {
        return true;
      }
      if (_top_peak_mz[i] < other._top_peak_mz[j]) {
        ++i;
      } else {
        ++j;
      }
    }
    return false;
  }
};

inline ostream& operator<<(ostream& os, const Spectrum& spectrum) {
  os << "is_clustered: " << spectrum._is_clustered << endl;
  os << "is_consensus: " << spectrum._is_consensus << endl;
  os << "precursor_mz: " << spectrum._precursor_mz << endl;
  os << "charge: " << spectrum._charge << endl;
  os << "peptide: " << spectrum._peptide << endl;
  os << "title: " << spectrum._title << endl;
  os << "component titles: " << spectrum._component_titles << endl;
  os << "raw peaks: " << endl;
  for (const auto& peak : spectrum._raw_peaks) {
    os << peak << endl;
  }
  os << endl;
  os << "filtered peaks: " << endl;
  for (const auto& peak : spectrum._filtered_peaks) {
    os << peak << endl;
  }
  os << endl;
  os << "top peaks' mz:" << endl;
  for (auto mz : spectrum._top_peak_mz) {
    os << mz << endl;
  }
  os << endl;
  os << "embeded_peaks: " << endl;
  for (const auto& peak : spectrum._embeded_peaks) {
    os << peak << endl;
  }
  os << endl;
  return os;
}

}  // namespace Core
#endif
