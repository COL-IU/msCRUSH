#ifndef __DISTANCE_H__
#define __DISTANCE_H__

#include <cmath>
#include <utility>
#include <unordered_map>
#include <vector>

#include "params.h"

using std::pair;
using std::unordered_map;
using std::vector;

namespace Core {
class Distance {
 public:
  static float euclidean(const EmbededPeaks& lhs,
      const EmbededPeaks& rhs);
  static float cosine(const EmbededPeaks& lhs,
      const EmbededPeaks& rhs);
  static float cosine(const Peaks& lhs, const Peaks& rhs, float precision);
  static float dot_product(const EmbededPeaks& lhs,
      const EmbededPeaks& rhs);

};

}  // namespace Core
#endif
