#include "io.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

/* Populates global famdist_map with root family distribution read from famdist_file_path */
void read_famdist(string famdist_file_path) {

  ifstream famdist_file(famdist_file_path.c_str()); // the constructor for ifstream takes const char*, not string, so we need to use c_str()
  string line;
  while (getline(famdist_file, line)) {
    cout << line << endl;
  }
}
