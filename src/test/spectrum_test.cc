#include <iostream>

#include "../class/spectrum.h"

using namespace Core;
using namespace std;

/*
ostream& operator<<(ostream& os, const Spectrum& spectrum) {
  os << "precursor_mz: " << spectrum._precursor_mz << endl;
  os << "charge: " << spectrum._charge << endl;
  os << "peptide: " << spectrum._peptide << endl;
  os << "title: " << spectrum._title << endl;
  os << "raw peaks: " << endl;
  for (const auto& peak : spectrum._raw_peaks) {
      os << peak._mz << " : " << peak._intensity << ", count: " << peak._count << endl;
  }
  os << endl;
  os << "filtered peaks: " << endl;
  for (const auto& peak : spectrum._filtered_peaks) {
      os << peak._mz << " : " << peak._intensity << endl; 
  }
  os << endl;
  os << "embeded_peaks: ";
  for (const auto& peak : spectrum._embeded_peaks) {
      os << peak._idx << " : " << peak._intensity << ", with count: " << peak._count << endl;
  }
  return os;
}*/

void TEST_Basics() {
  Spectrum spectrum;

  cout << "TEST uninitialized spectrum." << endl;
  cout << spectrum << endl;
  cout << endl;

  cout << "TEST initialized spectrum." << endl;
  spectrum._charge = 2;
  spectrum._precursor_mz = 100.;
  spectrum._peptide = "AGCT";
  spectrum._title = "AGCT TITLE";
  spectrum._raw_peaks.push_back(Peak(1.5, 1.5));
  spectrum._embeded_peaks.push_back(EmbededPeak(1, 1));
  cout << spectrum << endl;
}

void TEST_ShareTopPeaks_1() {
  Spectrum s1, s2;
  const float epsilon = 1;
  s1._top_peak_mz =
  {
    204.027, 
    204.951, 
    205.066, 
    211.019, 
    216, 
    218.006, 
    221.039, 
    222.962, 
    229.029, 
    230.024, 
  };
  s2._top_peak_mz = 
  {
  };
  cout << s1.shareTopPeaks(s2, epsilon) << endl;
}

void TEST_ShareTopPeaks_2() {
  Spectrum s1, s2;
  const float epsilon = 1;
  s2._top_peak_mz =
  {
    204.027, 
    204.951, 
    205.066, 
    211.019, 
    216, 
    218.006, 
    221.039, 
    222.962, 
    229.029, 
    230.024, 
  };
  s1._top_peak_mz = 
  {
  };
  cout << s1.shareTopPeaks(s2, epsilon) << endl;
}

void TEST_ShareTopPeaks_3() {
  Spectrum s1, s2;
  const float epsilon = 0.1;
  s2._top_peak_mz =
  {
    204.027, 
    204.951, 
    205.066, 
    211.019, 
    216, 
    218.006, 
    221.039, 
    222.962, 
    229.029, 
    230.024, 
  };
  s1._top_peak_mz = 
  {
      183.314,
      183.516,
      200.499,
      206.436,
      214.35 ,
      218.036,
      221.059,
  };
  cout << s1.shareTopPeaks(s2, epsilon) << endl;
}
int main() {
  TEST_Basics();
  TEST_ShareTopPeaks_3();
  return 0;
}
