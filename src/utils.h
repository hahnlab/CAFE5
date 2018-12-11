#ifndef utils_h
#define utils_h

#include <regex>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

class clade;

class newick_parser {

private:
    std::regex tokenizer;
    int lp_count, rp_count;

 /* methods */
    clade *new_clade(clade *p_parent);

public:
    bool parse_to_lambdas;	// flag for reading lambda tree
    std::string newick_string;

    /* methods */
    newick_parser(bool parse_lambdas) : tokenizer("\\(|\\)|[^\\s\\(\\)\\:\\;\\,]+|\\:[+-]?[0-9]*\\.?[0-9]+([eE][+-]?[0-9]+)?|\\,|\\;"),
        lp_count(0), rp_count(0), parse_to_lambdas(parse_lambdas) {} // constructor
    clade *parse_newick();
};

// these functions are intended to work with maps (key, value pairs)
template <typename T, typename U>
bool max_key(const std::pair<T, U> & p1, const std::pair<T, U> & p2) {
    return p1.first < p2.first;
}

template <typename T, typename U>
bool max_value(const std::pair<T, U> & p1, const std::pair<T, U> & p2) {
    return p1.second < p2.second;
}

//! Split string into vector of strings given delimiter
std::vector<std::string> tokenize_str(std::string some_string, char some_delim);
#endif
