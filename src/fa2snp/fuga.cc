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

void Fuga::Main(const string& f_nom) {
  ReadFaLst(f_nom);
  Print();
}

void Fuga::ReadFaLst(const string& f_nom) {
  ifstream f(f_nom);
  string l = "";
  std::vector<std::string> l_sp;
  getline(f, l);
  split(l_sp, l, "\t");
  SetCntgInfo(l_sp.at(1));
  f.seekg(0);
  while (getline(f, l)) {
    split(l_sp, l, "\t");
    iso.push_back(l_sp.at(0));
    ReadFa(l_sp.at(0), l_sp.at(1));
  }
}

void Fuga::SetCntgInfo(const string& f_nom) {
  ifstream f(f_nom);
  string l = "";
  std::vector<std::string> l_sp;
  string cntg_nom("");
  while (getline(f, l)) {
    if (l.size() == 0) {
      continue;
    } else if (l.at(0) == '>') {
      split(l_sp, l, " ");
      cntg_nom = l_sp.at(0).substr(1);
      cntg_info.insert(make_pair(cntg_nom, 0));
    } else if (l.at(0) != '>') {
      cntg_info.at(cntg_nom) += l.size();
    }
  }
}

void Fuga::ReadFa(const string& iso_nom, const string& f_nom) {
  cntg_seq.insert(make_pair(iso_nom, new unordered_map<string, string>()));
  ifstream f(f_nom);
  string l = "";
  std::vector<std::string> l_sp;
  string cntg_nom("");
  while (getline(f, l)) {
    if (l.size() == 0) {
      continue;
    } else if (l.at(0) == '>') {
      split(l_sp, l, " ");
      cntg_nom = l_sp.at(0).substr(1);
      (*cntg_seq.at(iso_nom)).insert(make_pair(cntg_nom, ""));
    } else {
      (*cntg_seq.at(iso_nom)).at(cntg_nom) += l;
    }
  }
}

void Fuga::Print() {
  cout << "contig" << '\t' << "pos";
  for (auto it = iso.begin(); it != iso.end(); ++it) {
    cout << '\t' << *it;
  }
  cout << '\n';
  for (auto it1 = cntg_info.begin(); it1 != cntg_info.end(); ++it1) {
    string cntg_nom = it1->first;
    for (ul pos = 0; pos < it1->second; ++pos) {
      unordered_set<char> nuc;
      for (auto it2 = cntg_seq.begin(); it2 != cntg_seq.end(); ++it2) {
        char _nuc((*it2->second).at(cntg_nom).at(pos));
        if (kValidNuc.find(_nuc) != string::npos) {
          nuc.insert(_nuc);
        }
      }
      if (nuc.size() >= 2) {
        cout << cntg_nom << '\t' << pos + 1;
        for (auto it2 = iso.begin(); it2 != iso.end(); ++it2) {
          cout << '\t' << (*cntg_seq.at(*it2)).at(cntg_nom).at(pos);
        }
        cout << '\n';
      }
    }
  }
  for (auto it = cntg_seq.begin(); it != cntg_seq.end(); ++it) {
    delete it->second;
  }
}
