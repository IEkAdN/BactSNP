#ifndef HOGE_H_
#define HOGE_H_

#include <algorithm>
#include <cctype>
#include <cfloat>
#include <climits>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::accumulate;
using std::cerr;
using std::cin;
using std::cout;
using std::greater;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::multimap;
using std::ofstream;
using std::pair;
using std::set;
using std::stod;
using std::stoi;
using std::stol;
using std::stoll;
using std::stoul;
using std::stoull;
using std::string;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

typedef int32_t i32;
typedef int64_t i64;
typedef uint32_t u32;
typedef uint64_t u64;

void split(vector<string>& s_sp, const string& s, const string& delim,
           bool token_compress = false);

void split_2(vector<string>& s_sp, const string& s, const string& delim,
           bool token_compress = false);

#endif  // HOGE_H_
