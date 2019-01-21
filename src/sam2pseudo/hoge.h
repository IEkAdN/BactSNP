/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#ifndef HOGE_H_
#define HOGE_H_

#include <algorithm>
#include <cfloat>
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
using std::stol;
using std::stoll;
using std::stoul;
using std::stoull;
using std::string;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::vector;

typedef unsigned long long ull;
typedef unsigned long ul;

void split(vector<string>& s_sp, const string& s, const string& delim,
           bool token_compress = false);

#endif  // HOGE_H_
