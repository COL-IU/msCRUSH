#ifndef __CORE_H__
#define __CORE_H__

#include <algorithm>
#include <random>
#include <string>
#include <vector>

#include "params.h"

using namespace std;

namespace Core {
typedef vector<float> hashFunction;
typedef vector<hashFunction> hashTable;
class LSH {
 public:
  static hashFunction generateNormalHashFunc(int size);
  static hashFunction generateUniformHashFunc(int size, float min, float max);
  static hashTable generateHashTable(int htable_size, int hfunction_size); 
  
  // Different LSHs.
  static int random_projection(const EmbededPeaks& embeded_peaks,
      const hashFunction& function, bool intensity_binarized = false);
  static int random_projection(const EmbededPeaks& embeded_peaks,
      const hashTable& table, bool intensity_binarized = false);
  static string p_stable(const EmbededPeaks& embeded_peaks,
      const hashFunction& function, float b, float r);
  static string p_stable(const EmbededPeaks& embeded_peaks,
      const hashTable& table, float b, float r);
};

}  // namespace Core
#endif
