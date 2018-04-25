#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
using namespace std;

int main (int argc, char *argv[]) {
  string file = argv[1];
  cout << "file name: " << file << endl;
  ifstream ifs(file);
  string line;
  while(getline(ifs, line));
}
