#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
using namespace std;

int main (int argc, char *argv[]) {
  string file = argv[1];
  int bin = stoi(argv[2]);
  cout << "file name: " << file << endl;
  ifstream is(file);
  is.seekg(0, is.end);
  cout << "total bytes: " << is.tellg() << endl;
  is.seekg(0, is.beg);

  string line;
  int spectrum_cnt = 0;
  vector<streampos> pos;
  streampos p;
  while(getline(is, line)) {
    if (string::npos != line.find("END")) {
      ++spectrum_cnt;
      if (spectrum_cnt != 0 && spectrum_cnt % bin == 0) {
        p = is.tellg();
        pos.push_back(p);
      }
    }
  }
  pos.push_back(is.tellg());

  cout << "token pos: " << endl;
  for (const auto& p : pos) {
    cout << (long long)p << endl;
    is.clear();
    is.seekg(p, is.beg);
    cout << is.tellg() << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
    getline(is, line);
    cout << "" << line << endl;
  }

  return 0;
}
