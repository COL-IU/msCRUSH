#include "io.h"
namespace Utility {

float IO::convert(char const* source, char ** endPtr ) {
  char* end;
  int left = strtol( source, &end, 10 );
  float results = left;
  if ( *end == '.' ) {
      char* start = end + 1;
      int right = strtol( start, &end, 10 );
      static double const fracMult[] 
          = { 0.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001 };
      results += right * fracMult[ end - start ];
  }
  if ( endPtr != nullptr ) {
      *endPtr = end;
  }
  return results;
}

void IO::ReadSpectraFromMGF(vector<Spectrum*>* indexed_spectra, 
    unordered_map<int, vector<int>>* map_spectra_by_charge,
    unordered_map<string, int>* map_ms_titles_to_index,
    int* spectra_size, string file_name, float scale, float min_mz,
    float max_mz, float precision, int topK, int bin_size,
    bool peak_normalized, bool remove_precursor, bool verbose) {

  // Index shall follow next MS/MS spectrum already stored in indexed_spectra.
  const int spectra_cnt_stored = (*indexed_spectra).size();
  int spectra_cnt = 0;

  Peaks raw_peaks (MAX_PEAK_SIZE);  // Pre-allocate memory.
  Peaks filtered_peaks;  // Filtered peaks using bins.
  int i_peak = 0;

  string line, tmp, title = "NA", peptide = "NA";
  int charge = 0;
  float mz = 0, intensity = 0, precursor_mz = 0;

  //TODO: 20ppm.  [abs(threotical - observed) / threotical] * 10^6.
  const float ppm = 20;
  float precursor_peak_tol = 0;
  float precursor_no_water_peak_tol = 0; 
  float peak_mz_remove_without_water = -1;


  ifstream infile(file_name);
  if (!infile.is_open()) {
    cout << "Error! file not open!" << endl;
    exit(-1);
  }

  while (getline(infile, line)) {
    // Handle text.
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (0 == line.find(PEPMASS_MGF)) {
      precursor_mz = stof(line.substr(PEPMASS_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(CHARGE_MGF)) {
      charge = stoi(line.substr(CHARGE_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(TITLE_MGF)) {
      title = line.substr(TITLE_MGF.length() + 1);

      // Get rid of '^M', carriage-and-return created by Windows platform.
      if (title.back() == '\r') {
        title.pop_back();
      }
      continue;
    }
    if (line[0] < '0' || line[0] > '9') {
      continue;
    }
    
    peak_mz_remove_without_water = 
      precursor_mz - 18./(charge == 0 ? 1 : charge);

    precursor_peak_tol =  precursor_mz * ppm / 1000000;
    precursor_no_water_peak_tol = peak_mz_remove_without_water * ppm / 1000000;

    // Initialize stuff. 
    // raw_peaks.clear();
    i_peak = 0;
    filtered_peaks.clear();
    
    // Read peaks.
    do {
      if (string::npos != line.find(END_IONS_MGF)) {
        break;
      }
      if (line.empty()) {
        continue;
      }

      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      char* end;
      mz = convert(line.c_str(), &end);
      intensity = convert(end, &end);

      // Take the natural logarithm of intensity, mitigating dominant peaks.
      if (peak_normalized) {
        intensity = intensity == 0 ? 0 : log(intensity);
      }

      // Filter out peaks out of mz range.
      if (mz < min_mz || mz > max_mz || intensity <= 0) {
        continue;
      }

      // Filter out peaks around precursor_mz and precursor_mz- H2O/charge.
      if (remove_precursor && (
          fabs(mz - precursor_mz) < precursor_peak_tol ||
          fabs(mz - peak_mz_remove_without_water) <
          precursor_no_water_peak_tol)) {
        continue;
      }

      // raw_peaks.push_back(Peak(mz, intensity));
      raw_peaks[i_peak++] = Peak(mz, intensity);

    } while(getline(infile, line)); 

    BinTopKPeak(&filtered_peaks, raw_peaks, i_peak, topK, bin_size);
    RemoveAdjacentPeaks(&filtered_peaks, precision);
    Normalize(&filtered_peaks, scale);
    EmbededPeaks embeded_peaks;
    Embed(&embeded_peaks, filtered_peaks, min_mz, precision, scale);

    vector<float> top_peak_mz = SelectTopPeakMZ(filtered_peaks);

    // For memory efficient. Discard these 2 lines if memory is sufficient.
    // If descarded, have to clear these vectors before @51L-52L.
    // raw_peaks.clear();
    // filtered_peaks.clear();

    Spectrum* spectrum = new Spectrum();
    SetSpectrum(spectrum, false, false, charge, 1, embeded_peaks, Peaks(),
        filtered_peaks, peptide, precursor_mz, title, title, top_peak_mz);

    (*indexed_spectra).push_back(spectrum);
    
    (*map_spectra_by_charge)[charge].push_back(
        spectra_cnt_stored + spectra_cnt);
    (*map_ms_titles_to_index)[title] = spectra_cnt_stored + spectra_cnt;
    
    ++spectra_cnt;
    if (verbose && 0 == spectra_cnt % 100000) {
      cout << "read #spectra: " << spectra_cnt << endl;
    }

    // Reset charge, peptide, title, precursor_mz, intensity.
    charge = 0;
    precursor_mz = 0;
    intensity = 0;
    peptide = "NA";
    title = "NA";
  }
  //cout << "read #spectra: " << spectra_cnt << endl;

  (*spectra_size) = spectra_cnt;
  //cout << "read #spectra: " << spectra_cnt << ", from: " << file_name << endl;
  infile.close();
}

vector<float> IO::SelectTopPeakMZ(const Peaks& _peaks, int topK) {
  Peaks peaks(_peaks);
  auto cmp = [](const Peak& p1,const Peak& p2)
  {return p1._intensity > p2._intensity;};

  topK = min(topK, (int)peaks.size());
  nth_element(peaks.begin(), peaks.begin() + topK - 1, peaks.end(), cmp);
  vector<float> top_peak_mz(topK);
  for (int i = 0; i < topK; ++i) {
    top_peak_mz[i] = peaks[i]._mz;
  }
  // Sort according to mz ascendingly.
  sort(top_peak_mz.begin(), top_peak_mz.end());
  return top_peak_mz;
}

void IO::SetSpectrum(Spectrum* spectrum,
    bool is_clustered, bool is_consensus, int charge,
    int count, const EmbededPeaks& embeded_peaks, const Peaks& raw_peaks,
    const Peaks& filtered_peaks, string peptide, float pre_mz,
    string title, string component_titles, const vector<float>& top_peak_mz) {
  (*spectrum)._charge = charge;
  (*spectrum)._component_titles = component_titles;
  (*spectrum)._count = count;
  (*spectrum)._embeded_peaks = embeded_peaks;
  (*spectrum)._filtered_peaks = filtered_peaks;
  (*spectrum)._is_clustered = is_clustered;
  (*spectrum)._is_consensus = is_consensus;
  (*spectrum)._peptide = peptide;
  (*spectrum)._precursor_mz = pre_mz;
  (*spectrum)._raw_peaks = raw_peaks;
  (*spectrum)._title = title;
  (*spectrum)._top_peak_mz = top_peak_mz;
}

void IO::RemoveAdjacentPeaks(Peaks* peaks, float mz_tolerance) {
  if ((*peaks).empty()) {
    return;
  }
  int i_read, i_write;
  for (i_write = 0, i_read = 1; i_read < (int)(*peaks).size(); ++i_read) {
    auto r_spectrum = (*peaks)[i_read];
    auto& w_spectrum = (*peaks)[i_write]; 
    if (fabs(r_spectrum._mz - w_spectrum._mz) <= mz_tolerance) {
      if (r_spectrum._intensity > w_spectrum._intensity) {
        w_spectrum = r_spectrum;
      }
    } else {
      (*peaks)[++i_write] = r_spectrum;
    }
  }
  (*peaks).resize(i_write + 1);
}

void IO::Normalize(Peaks* peaks, float scale) {
  float max_intensity = 0;
  for (const auto& peak : *peaks) {
    max_intensity = max(max_intensity, peak._intensity);
  }
  for (auto& peak : *peaks) {
    peak._intensity = peak._intensity / max_intensity * scale;
  }
}

void IO::Embed(EmbededPeaks* embeded_peaks, const Peaks& peaks,
    float min_mz, float precision, float scale) {
  (*embeded_peaks).clear();
  for (const auto& peak : peaks) {
    float mz = peak._mz;
    float intensity = peak._intensity;
    int idx = round((mz - min_mz) / precision);
    (*embeded_peaks).push_back(EmbededPeak(idx, intensity, peak._count));
  }
}

void IO::MergeTwoPeaks (const Peaks& peaks1, const Peaks& peaks2, Peaks* peaks) {
  
  int i = 0, j = 0;
  int i_len = peaks1.size(), j_len = peaks2.size();
  Peaks local_peaks;
  local_peaks.reserve(i_len + j_len);

  while (i < i_len && j < j_len) {
    if (peaks1[i]._mz == peaks2[j]._mz) {
      //TODO: Need once or twice?
      local_peaks.push_back(peaks1[i]);
      local_peaks.push_back(peaks2[j]);
      ++i;
      ++j;
    } else if (peaks1[i]._mz < peaks2[j]._mz) {
      local_peaks.push_back(peaks1[i]);
      ++i;
    } else {
      local_peaks.push_back(peaks2[j]);
      ++j;
    }
  }

  while (i < i_len) {
    local_peaks.push_back(peaks1[i]);
    ++i;
  }

  while (j < j_len) {
    local_peaks.push_back(peaks2[j]);
    ++j;
  }

  //*peaks = move(local_peaks);
  swap(*peaks, local_peaks);
}

// Merge spectra to build a consensus spectrum.
void IO::SetConsensus(Spectrum* consensus,
    const Spectrum& s1, const Spectrum& s2,
    float precision, int topK, float bin_size, float min_mz, float scale, 
    string title, string component_titles) {

  Peaks all_peaks, peaks1, peaks2;
  peaks1.reserve(s1._filtered_peaks.size());
  peaks2.reserve(s1._filtered_peaks.size());

  Peak tmp;
  int count1 = s1._count;
  for (const auto& peak : s1._filtered_peaks) {
    tmp = peak;
    tmp._intensity *= count1;
    peaks1.push_back(tmp);
  }

  int count2= s2._count;
  for (const auto& peak : s2._filtered_peaks) {
    tmp = peak;
    tmp._intensity *= count2;
    peaks2.push_back(tmp);
  }

  MergeTwoPeaks(peaks1, peaks2, &all_peaks);

  int total_count = count1 + count2;
  float ave_precursor_mz =
    (s1._precursor_mz * count1 + s2._precursor_mz * count2 ) / total_count;

  int size = s1._filtered_peaks.size() + s2._filtered_peaks.size();

  //TODO: Better options? Can change to while loop for more iterations.
  Peaks new_peaks(size);
  float threshold = precision;
  new_peaks[0] = all_peaks[0];
  int i_last = 0;
  for (int i = 1; i < all_peaks.size(); ++i) {
    if (all_peaks[i]._mz - new_peaks[i_last]._mz <= threshold) {
      auto& last_peak = new_peaks[i_last];
      const auto& to_merge_peak = all_peaks[i];
      float sum_inten = last_peak._intensity + to_merge_peak._intensity;

      // Set summed inten, weighted mz, count, i.e(spectrum that has it peak)
      last_peak._mz = last_peak._mz * last_peak._intensity / sum_inten +
        to_merge_peak._mz * to_merge_peak._intensity / sum_inten;
      last_peak._intensity = sum_inten;

    } else {
      // new_peaks.push_back(all_peaks[i]); 
      new_peaks[++i_last] = all_peaks[i];
    }
  }
  new_peaks.resize(++i_last);

  Peaks filtered_peaks;
  BinTopKPeak(&filtered_peaks, new_peaks, i_last, topK, bin_size);
  RemoveAdjacentPeaks(&filtered_peaks, precision);
  Normalize(&filtered_peaks, scale);
  EmbededPeaks embeded_peaks;
  Embed(&embeded_peaks, filtered_peaks, min_mz, precision, scale);

  vector<float> top_peak_mz = SelectTopPeakMZ(filtered_peaks);

  // For memory efficient. Discard these 2 lines if memory is sufficient.
  new_peaks.clear();
  // filtered_peaks.clear();

  int charge = (0 == s1._charge ? s2._charge: s1._charge);
  SetSpectrum(consensus, false, true, charge, total_count, embeded_peaks,
      new_peaks, filtered_peaks, "CONSENSUS_SPECTRUM_PEPTIDE",
      ave_precursor_mz, title, component_titles, top_peak_mz);

}


// Adapt the peak intensities in consensus peak using the following formula:
// I = I * (0.95 + 0.05 * (1 + pi)^5)
void IO::AdaptPeakIntensities(Peaks* peaks, int nSpectra) {
  for (auto& peak : *peaks) {
    float pi = (float)peak._count / nSpectra;
    //TODO: Plan to deprecate it.
    // assert(pi <= 1);
    if (pi > 1) {
      cout << peak << endl;
      cout << "nSpectra: " << nSpectra << endl;
      cout << "Pi greater than 1!" << endl;
      pi = 1.;
    }
    float new_intensity = peak._intensity * (0.95 + 0.05 * pow(1 + pi, 5));
    peak._intensity = new_intensity;
  }
}

void IO::BinTopKPeak(Peaks* top_peaks, const Peaks& peaks, int peaks_size,
    int topK, float bin_size) {
  // Assume that peaks were already sorted ascendingly w.r.t mz.
  float max_mz = peaks[peaks_size - 1]._mz;
  // Sort peaks descendingly w.r.t peak's intensity.
  auto cmp =
    [](const Peak& p1, const Peak& p2){return p1._intensity > p2._intensity;};

  *top_peaks = Peaks((round(max_mz / bin_size) + 1 ) * topK);
  Peaks current_peaks(peaks_size);
  int i_write = 0;

  int last_index = 0;
  float current_max_mz = 0;
  for (float min_mz = 0; min_mz < max_mz; min_mz += bin_size) {
    current_max_mz = min_mz + bin_size;
    int i_current_peak = 0;
    for (; last_index < peaks_size; ++last_index) {
      const auto& peak = peaks.at(last_index);
      if (peak._mz > current_max_mz) {
        break;
      }
      current_peaks[i_current_peak++] = peak;
    }

    int ele_num = min(i_current_peak, topK);
    //TODO: Should be .begin() + ele_num - 1.
    nth_element(current_peaks.begin(), current_peaks.begin() + ele_num - 1,
        current_peaks.begin() + i_current_peak, cmp);   
    for (int i = 0; i < ele_num; ++i) {
      (*top_peaks)[i_write++] = current_peaks[i];
    }
  }
  (*top_peaks).resize(i_write);
  sort((*top_peaks).begin(), (*top_peaks).end());
}

}  // namespace Utility
