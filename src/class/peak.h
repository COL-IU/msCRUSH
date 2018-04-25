#ifndef __CORE_PEAK_H__
#define __CORE_PEAK_H__

#include <iostream>

using std::endl;
using std::ostream;

namespace Core {
class Peak {
 public:
  float _mz;
  float _intensity;
  int _count;

  Peak(): _mz(-1), _intensity(-1), _count(-1){}
  ~Peak() {}

  Peak(const Peak& p) {
    _count = p._count;
    _intensity = p._intensity;
    _mz = p._mz;
  }
  Peak(float mz, float intensity, int count = 1) {
    _count = count;
    _intensity = intensity;
    _mz = mz;
  }

  Peak& operator=(const Peak& p) {
    _count = p._count;
    _intensity = p._intensity;
    _mz = p._mz;
    return *this;
  }

  bool operator<(const Peak& p) {
    return _mz < p._mz;
  }

  friend ostream& operator<<(ostream& os, const Peak& peak);
};

inline ostream& operator<<(ostream& os, const Peak& peak) {
  os << "[Peak] mz: " << peak._mz;
  os << ", intensity: " << peak._intensity;
  os << ", count: " << peak._count;
  return os;
}

class EmbededPeak {
 public:
  float _intensity;
  int _count;
  int _idx;

  EmbededPeak(): _idx(-1), _intensity(-1), _count(-1){}
  ~EmbededPeak() {}

  EmbededPeak(const EmbededPeak& p) {
    _count = p._count;
    _idx = p._idx;
    _intensity = p._intensity;
  }
  EmbededPeak(int idx, float intensity, int count = 1) {
    _count = count;
    _idx = idx;
    _intensity = intensity;
  }

  EmbededPeak& operator=(const EmbededPeak& p) {
    _count = p._count;
    _idx = p._idx;
    _intensity = p._intensity;
    return *this;
  }

  friend ostream& operator<<(ostream& os, const EmbededPeak& peak);
};

inline ostream& operator<<(ostream& os, const EmbededPeak& peak) {
  os << "[EmbededPeak] idx: " << peak._idx;
  os << ", intensity: " << peak._intensity;
  os << ", count: " << peak._count;
  return os;
}

}  // namespace Core
#endif
