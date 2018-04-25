#include "distance.h"
namespace Core{

//TODO: Implement a highly efficient method like 'cosine'.
float Distance::euclidean(const EmbededPeaks& lhs,
  const EmbededPeaks& rhs) {
  float distance = 0;
  unordered_map<int, double> rhs_map;
  unordered_map<int, bool> rhs_visited; //  Idx is visited or not

  for (const auto & peak: rhs) {
    rhs_map[peak._idx] = peak._intensity;
    rhs_visited[peak._idx] = false;
  }
  
  for (const auto & peak : lhs) {
    if(rhs_map.find(peak._idx) != rhs_map.end()) {
      distance += pow(peak._intensity - rhs_map[peak._idx], 2);
      rhs_visited[peak._idx] = true;
    } else {
      distance += pow(peak._intensity, 2);
    }
  }
      
  for (const auto & peak : rhs) {
    if (!rhs_visited[peak._idx]) {
      distance += pow(peak._intensity, 2);
    }
  }

  return sqrt(distance);
}

float Distance::cosine(const Peaks& lhs, const Peaks& rhs, float precision) {
  float similarity = 0;
  int i = 0, j = 0;
  float i_mz, j_mz;
  float magnitude_lhs = 0, magnitude_rhs = 0;
  float i_intensity, j_intensity;
  while (i < lhs.size() && j < rhs.size()) {
    i_mz = lhs[i]._mz;
    j_mz = rhs[j]._mz;

    i_intensity = lhs[i]._intensity;
    j_intensity = rhs[j]._intensity;

    if (fabs(i_mz - j_mz) <= precision) {
      similarity += i_intensity * j_intensity;
      magnitude_lhs += i_intensity * i_intensity;
      magnitude_rhs += j_intensity * j_intensity;
      ++i;
      ++j;
    } else if (i_mz < j_mz) {
      magnitude_lhs += i_intensity * i_intensity;
      ++i;
    } else {
      magnitude_rhs += j_intensity * j_intensity;
      ++j;
    }
  }
  while (i < lhs.size()) {
    i_intensity = lhs[i++]._intensity;
    magnitude_lhs += i_intensity * i_intensity;
  }

  while (j < rhs.size()) {
    j_intensity = rhs[j++]._intensity;
    magnitude_rhs += j_intensity * j_intensity;
  }
  
  similarity /= sqrt(magnitude_lhs) * sqrt(magnitude_rhs);
  return 1 - similarity;
}
  
float Distance::cosine(const EmbededPeaks& lhs,
        const EmbededPeaks& rhs) {
  float similarity = 0;
  int i = 0, j = 0, i_idx, j_idx;
  float magnitude_lhs = 0, magnitude_rhs = 0;
  float i_intensity, j_intensity;
  while (i < lhs.size() && j < rhs.size()) {
    i_idx = lhs[i]._idx;
    j_idx = rhs[j]._idx;
    i_intensity = lhs[i]._intensity;
    j_intensity = rhs[j]._intensity;
    if (i_idx == j_idx) {
      similarity += i_intensity * j_intensity;
      magnitude_lhs += i_intensity * i_intensity;
      magnitude_rhs += j_intensity * j_intensity;
      ++i;
      ++j;
    } else if (i_idx < j_idx) {
      magnitude_lhs += i_intensity * i_intensity;
      ++i;
    } else {
      magnitude_rhs += j_intensity * j_intensity;
      ++j;
    }
  }

  while (i < lhs.size()) {
    i_intensity = lhs[i++]._intensity;
    magnitude_lhs += i_intensity * i_intensity;
  }

  while (j < rhs.size()) {
    j_intensity = rhs[j++]._intensity;
    magnitude_rhs += j_intensity * j_intensity;
  }
  
  similarity /= sqrt(magnitude_lhs) * sqrt(magnitude_rhs);
  return 1 - similarity;
}

//TODO: refactor dot_product from github repository. 
//https://github.iu.edu/wang558/lsh_clustering/blob/master/src/nearest_neightbor.cpp
float Distance::dot_product(const EmbededPeaks& lhs,
        const EmbededPeaks& rhs) {
  return 0;
}

}  //namespace Core
