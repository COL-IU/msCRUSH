#include <cassert>
#include <iostream>
using namespace std;

#include "../class/distance.h"
using namespace Core;

const float eps = 0.00001;

void TEST_COSINE1() {
  EmbededPeaks lhs {EmbededPeak(1, 2.), EmbededPeak(2, 3.), EmbededPeak(4, -1.)};
  EmbededPeaks rhs {EmbededPeak(1, 2.), EmbededPeak(3, 3.), EmbededPeak(4, 0.5)};
  float dist = Distance::cosine(lhs, rhs);
  assert(abs(dist - 0.743022) < eps);
  cout << "TEST_COSINE1 success!" << endl;
}

void TEST_COSINE2() {
  EmbededPeaks lhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  EmbededPeaks rhs {EmbededPeak(1, 8.), EmbededPeak(2, 1.), EmbededPeak(3, -3), EmbededPeak(4, 4)};
  float dist = Distance::cosine(lhs, rhs);
  assert(abs(dist - 1.138013) < eps);
  cout << "TEST_COSINE2 success!" << endl;
}

void TEST_COSINE3() {
  EmbededPeaks lhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  EmbededPeaks rhs {EmbededPeak(1, 1)};
  float dist = Distance::cosine(lhs, rhs);
  assert(abs(dist - 1) < eps);
  cout << "TEST_COSINE3 success!" << endl;
}

void TEST_COSINE4() {
  EmbededPeaks rhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  EmbededPeaks lhs {EmbededPeak(1, 8.), EmbededPeak(2, 1.), EmbededPeak(3, -3), EmbededPeak(4, 4)};
  float dist = Distance::cosine(lhs, rhs);
  assert(abs(dist - 1.138013) < eps);
  cout << "TEST_COSINE4 success!" << endl;
}

void TEST_EUCLIDEAN1() {
  vector<EmbededPeak> lhs {EmbededPeak(1, 2.), EmbededPeak(2, 3.), EmbededPeak(4, -1.)};
  vector<EmbededPeak> rhs {EmbededPeak(1, 2.), EmbededPeak(3, 3.), EmbededPeak(4, 0.5)};
  float dist = Distance::euclidean(lhs, rhs);
  assert(abs(dist - 4.5) < eps);
  cout << "TEST_EUCLIDEAN1 success!" << endl;
}

void TEST_EUCLIDEAN2() {
  vector<EmbededPeak> lhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  vector<EmbededPeak> rhs {EmbededPeak(1, 8.), EmbededPeak(2, 1.), EmbededPeak(3, -3), EmbededPeak(4, 4)};
  float dist = Distance::euclidean(lhs, rhs);
  assert(abs(dist - 10.062306) < eps);
  cout << "TEST_EUCLIDEAN2 success!" << endl;
}

void TEST_EUCLIDEAN3() {
  vector<EmbededPeak> lhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  vector<EmbededPeak> rhs {EmbededPeak(1, 1)};
  float dist = Distance::euclidean(lhs, rhs);
  assert(abs(dist - 2.5) < eps);
  cout << "TEST_EUCLIDEAN3 success!" << endl;
}

void TEST_EUCLIDEAN4() {
  vector<EmbededPeak> rhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  vector<EmbededPeak> lhs {EmbededPeak(1, 8.), EmbededPeak(2, 1.), EmbededPeak(3, -3), EmbededPeak(4, 4)};
  float dist = Distance::euclidean(lhs, rhs);
  assert(abs(dist - 10.062306) < eps);
  cout << "TEST_EUCLIDEAN4 success!" << endl;
}

void TEST_COSINE_IDX() {
  // Test cosine distance.
  TEST_COSINE1();
  TEST_COSINE2();
  TEST_COSINE3();
  TEST_COSINE4();

  TEST_EUCLIDEAN1();
  TEST_EUCLIDEAN2();
  TEST_EUCLIDEAN3();
  TEST_EUCLIDEAN4();
}

//TODO(wang558): TEST cosine precision
void TEST_COSINE_PRECISION() {
  return;
}

int main() {
  // TEST_COSINE_IDX();

  return 0;
}
