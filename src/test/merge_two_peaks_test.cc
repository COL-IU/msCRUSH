#include <iostream>

#include "../utility/io.h"
using namespace std;
using namespace Core;
using namespace Utility;


int main (int argc, char *argv[]) {

  Peaks peaks1 {Peak(1, 10,2), Peak(2.5, 10,2), Peak(2.8, 10,2), Peak(3.5, 10,2)};
  Peaks peaks2 {Peak(1.2, 10,2), Peak(1.5, 10,2), Peak(2.7, 10,2), Peak(3.5, 10,2)};

  Peaks peaks;
  IO::MergeTwoPeaks(peaks1, peaks2, &peaks);

  for (const auto& peak : peaks)
    cout << peak << endl;
  return 0;
}
