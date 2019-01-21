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

#include "fuga.h"

void Fuga::Main(string ref_nom, string region_nom) {
  ReadFa(ref_nom);
  ReadRegion(region_nom);
  Print();
}

void Fuga::ReadFa(string f_nom) {
  ifstream f(f_nom);
  string l;
  string header;
  while (getline(f, l)) {
    if (! l.empty()) {
      if (l.at(0) == '>') {
        header = l.substr(1, l.find(" ") - 1);
        fa_header_.push_back(header);
        fa_[header] = "";
      } else {
        fa_.at(header) += l;
      }
    }
  }
  for (auto it = fa_header_.begin(); it != fa_header_.end(); ++it) {
    region_[*it].resize(fa_.at(*it).size(), false);
  }
}

void Fuga::ReadRegion(string f_nom) {
  ifstream f(f_nom);
  string l;
  while (getline(f, l)) {
    vector<string> l_sp;
    split(l_sp, l, "\t");
    string cntg_nom(l_sp.at(0));
    ul beg(stoul(l_sp.at(1)) - 1);
    ul end(stoul(l_sp.at(2)) - 1);
    for (ul i = beg; i <= end; ++i) {
      region_.at(cntg_nom).at(i) = true;
    }
  }
}

void Fuga::Print() {
  for (auto it = fa_header_.begin(); it != fa_header_.end(); ++it) {
    string cntg(*it);
    bool pre_bool(false);
    for (ul pos = 0; pos < region_.at(cntg).size(); ++pos) {
      bool cur_bool(region_.at(cntg).at(pos));
      if ((! pre_bool) && cur_bool) {
        cout << cntg << "\t" << pos + 1;
      } else if (pre_bool && (! cur_bool)) {
        cout << "\t" << pos << "\n";
      }
      pre_bool = cur_bool;
    }
    if (pre_bool) {
      cout << "\t" << region_.at(cntg).size() << "\n";
    }
  }
}
