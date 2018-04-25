#include <iostream>
#include <cassert>

#include "../class/core.h"

using namespace Core;
using namespace std;

void TEST_Normal() {
  int size = 1000, min = 0, max = 0;
  hashFunction h = LSH::generateNormalHashFunc(size);
  for (const auto& ele : h) {
    cout << ele << endl;
  }
}

void TEST_generateHashTable() {
  int htable_size = 2, hfunction_size = 10;
  hashTable table = LSH::generateHashTable(htable_size, hfunction_size);
  cout << table.size() << endl;
  int cnt = 0;
  for (auto& f : table) {
    cout << "Table: " << ++cnt  << endl;
    for (auto& ele : f) {
      cout << ele << "\t";
    }
    cout << endl;
  }
}

void TEST_RANDOMPROJECTION() {
  EmbededPeaks rhs {EmbededPeak(1, 0.), EmbededPeak(2, 1.), EmbededPeak(3, 2), EmbededPeak(4, 0.5)};
  EmbededPeaks lhs {EmbededPeak(1, 8.), EmbededPeak(2, 1.), EmbededPeak(3, -3), EmbededPeak(7, 4)};
  hashTable table = LSH::generateHashTable(2, 10);

  int cnt = 0;
  for (auto& f : table) {
    cout << "Table: " << ++cnt  << endl;
    for (auto& ele : f) {
      cout << ele << endl;
    }
    cout << endl;
  }

  int key = LSH::random_projection(lhs, table);
  cout << "key: " << key << endl;
  cout << "TEST_RANDOMPROJECTION success!" << endl;
}

int main () {
  // TEST_Normal();
  // TEST_generateHashTable();
  TEST_RANDOMPROJECTION();
  return 0;
}
