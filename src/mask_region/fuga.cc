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

void Fuga::Main(string fa_nom, string pos_nom) {
  ReadFa(fa_nom);
  ReadRegion(pos_nom);
  Print();
}

void Fuga::ReadFa(string f_nom) {
  ifstream f(f_nom);
  string l = "";
  vector<string> l_sp;
  string cntg_nom("");
  while (getline(f, l)) {
    if (l.size() == 0) {
      continue;
    } else if (l.at(0) == '>') {
      split(l_sp, l, " ");
      cntg_nom = l_sp.at(0).substr(1);
      cntg_nom_.push_back(cntg_nom);
      cntg_seq_.insert(make_pair(cntg_nom, ""));
    } else {
      cntg_seq_.at(cntg_nom) += l;
    }
  }
}

void Fuga::ReadRegion(string f_nom) {
  ifstream f(f_nom);
  string l = "";
  vector<string> l_sp;
  while (getline(f, l)) {
    split(l_sp, l, "\t");
    string cntg_nom(l_sp.at(0));
    ul beg(stoul(l_sp.at(1)) - 1);
    ul end(stoul(l_sp.at(2)) - 1);
    ul len(end - beg + 1);
    cntg_seq_.at(cntg_nom).replace(beg, len, len, '-');
  }
}

void Fuga::Print() {
  for (auto it = cntg_nom_.begin(); it != cntg_nom_.end(); ++it) {
    const string& nom(*it);
    const string& seq(cntg_seq_.at(*it));
    cout << '>' << nom << '\n';
    for (ul i = 0; i < seq.size(); ++i) {
      cout << seq.at(i);
      if (i % 80 == 79) {
        cout << '\n';
      }
    }
    if (seq.size() % 80 != 0) {
      cout << '\n';
    }
  }
}
