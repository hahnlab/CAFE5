#ifndef io_h
#define io_h

#include <string>
#include <map>

using namespace std;

extern struct option longopts[];

map<int, int> *read_rootdist(string famdist_file_path);

#endif
