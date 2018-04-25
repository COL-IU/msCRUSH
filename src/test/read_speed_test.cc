#include <chrono>
#include <algorithm>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>  
#include <unistd.h>
#include "../utility/params.h"
using namespace std;
using namespace Utility;

float convert(char const* source, char ** endPtr) {
  char* end;
  int left = strtol( source, &end, 10 );
  float results = left;
  if ( *end == '.' ) {
      char* start = end + 1;
      int right = strtol( start, &end, 10 );
      static double const fracMult[] 
          = { 0.0, 0.1, 0.01, 0.001, 0.0001, 0.00001 };
      results += right * fracMult[ end - start ];
  }
  if ( endPtr != nullptr ) {
      *endPtr = end;
  }
  return results;
}

//struct membuf
//    std::streambuf {
//    membuf(char* start, size_t size) {
//        this->setg(start, start, start + size);
//    }
//};

void test_read_plain_text(string file) {
  cout << file << endl;
  ifstream iss(file);
  if (!iss.is_open()) {
    cout << "file not open!" << endl;
    exit(-1);
  }

  string line;
  long long line_num = 0;

  auto start_time = chrono::high_resolution_clock::now();
  iss.clear();
  iss.seekg(0, iss.beg);
  while (getline(iss, line));
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
            end_time - start_time).count();
  cout << "read file takes secs: " << elapsed_read << endl;


  line_num = 0;
  start_time = chrono::high_resolution_clock::now();
  iss.clear();
  iss.seekg(0, iss.beg);
  while (getline(iss, line)) {
    ++line_num;
  }
  end_time = chrono::high_resolution_clock::now();
  elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
            end_time - start_time).count();
  cout << "read file & cnt line number takes secs: " << elapsed_read << endl;

  cout << "#line: " << line_num << endl;
  iss.close();
}

void test_read_spectrum(string file) {
  string line, tmp, title, peptide;
  int charge;
  float mz, intensity, precursor_mz;
  float min_mz = 200, max_mz = 2000;
  ifstream infile(file);
  if (!infile.is_open()) {
    cout << "file not open!" << endl;
    exit(-1);
  }

  auto start_time = chrono::high_resolution_clock::now();
  long long cnt = 0;
  while (getline(infile, line)) {
    // Handle text.
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (0 == line.find(PEPMASS_MGF)) {
      precursor_mz = stof(line.substr(PEPMASS_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(CHARGE_MGF)) {
      charge = stoi(line.substr(CHARGE_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(TITLE_MGF)) {
      title = line.substr(TITLE_MGF.length() + 1);
      continue;
    }
    if (!isdigit(line[0])) {
      continue;
    }
    
    // Read peaks.
    do {
      if (line.empty()) {
        continue;
      }
      if (line[0] == 'E') {
        break;
      }

      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      char* end;
      //mz = strtof(line.c_str(), &end);
      //intensity = strtof(end, &end);
      mz = convert(line.c_str(), &end);
      intensity = convert(end, &end);

      // Filter out peaks with intensity out of mz range.
      if (mz < min_mz || mz > max_mz) {
        continue;
      }
      // raw_peaks.push_back(Peak(mz, intensity));
    } while(getline(infile, line)); 
  }
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
            end_time - start_time).count();
  cout << "read file takes secs: " << elapsed_read << endl;
  infile.close(); 
}

void test_mmap(string file) {
  string line, tmp, title, peptide;
  int charge;
  float mz, intensity, precursor_mz;
  float min_mz = 200, max_mz = 2000;
  int fd = open(file.c_str(), O_RDONLY);
  size_t len = lseek(fd,0,SEEK_END);
  char *buf = (char *) mmap(NULL,len,PROT_READ,MAP_PRIVATE,fd,0);    
  cout << "#bytes read: " << len << endl;
  istringstream infile(buf);

  auto start_time = chrono::high_resolution_clock::now();

  while (getline(infile, line)) {
    // Handle text.
    if (line.empty() || line[0] == '#') {
      continue;
    }
    if (0 == line.find(PEPMASS_MGF)) {
      precursor_mz = stof(line.substr(PEPMASS_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(CHARGE_MGF)) {
      charge = stoi(line.substr(CHARGE_MGF.length() + 1));
      continue;
    }
    if (0 == line.find(TITLE_MGF)) {
      title = line.substr(TITLE_MGF.length() + 1);
      continue;
    }
    if (!isdigit(line[0])) {
      continue;
    }
    
    // Read peaks.
    do {
      if (line.empty()) {
        continue;
      }
      if (line[0] == 'E') {
        break;
      }

      // https://stackoverflow.com/questions/17465061/how-to-parse-space-separated-floats-in-c-quickly
      char* end;
      //mz = strtof(line.c_str(), &end);
      //intensity = strtof(end, &end);
      mz = convert(line.c_str(), &end);
      intensity = convert(end, &end);

      // Filter out peaks with intensity out of mz range.
      if (mz < min_mz || mz > max_mz) {
        continue;
      }
      // raw_peaks.push_back(Peak(mz, intensity));
    } while(getline(infile, line)); 
  }
  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
            end_time - start_time).count();
  cout << "read file takes secs: " << elapsed_read << endl;
  munmap(buf, len);
  close(fd);
}

int main (int argc, char *argv[]) {
  // test_read_plain_text(argv[1]);
  // test_read_spectrum(argv[1]);
  
  test_mmap(argv[1]);

  return 0;
}
