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

void Fuga::Main(string snp_nom, string ref_nom, string out_dir_nom) {
  ReadRef(ref_nom);
  ReadSnp(snp_nom);
  Print(out_dir_nom);
}

void Fuga::ReadRef(string f_nom) {
  ifstream f(f_nom);
  string l;
  string header;
  while (getline(f, l)) {
    if (! l.empty()) {
      if (l.at(0) == '>') {
        header = l.substr(1, l.find(" ") - 1);
        cntg_id_.push_back(header);
        ref_[header] = "";
      } else {
        ref_.at(header) += l;
      }
    }
  }
}

void Fuga::ReadSnp(string f_nom) {
  ifstream f(f_nom);
  string l("");
  getline(f, l);
  vector<string> l_sp;
  split(l_sp, l, "\t");
  for (unsigned i = 2; i < l_sp.size(); ++i) {
    const string& iso(l_sp.at(i));
    iso_.push_back(iso);
    for (auto cntg_id = cntg_id_.begin(); cntg_id != cntg_id_.end(); ++cntg_id) {
      cntg_seq_[iso][*cntg_id] = ref_.at(*cntg_id);
    }
  }
  while (getline(f, l)) {
    split(l_sp, l, "\t");
    if (! IsAmbiguous(l_sp)) {
      const string& cntg_id(l_sp.at(0));
      ul pos(stoul(l_sp.at(1)) - 1);
      for (unsigned i = 0; i < iso_.size(); ++i) {
        const string& iso(iso_.at(i));
        cntg_seq_.at(iso).at(cntg_id).replace(pos, 1, 1, l_sp.at(i + 2).at(0));
      }
    }
  }
}

bool Fuga::IsAmbiguous(const vector<string>& l_sp) {
  string acgt("ACGT");
  for (auto i = l_sp.begin() + 2; i < l_sp.end(); ++i) {
    if (acgt.find(*i) == string::npos) {
      return true;
    }
  }
  return false;
}

void Fuga::Print(string dir_nom) {
  for (auto iso = iso_.begin(); iso != iso_.end(); ++iso) {
    ofstream fa(dir_nom + "/" + *iso + ".fa");
    for (auto cntg_id = cntg_id_.begin(); cntg_id != cntg_id_.end(); ++cntg_id) {
      fa << ">" << *cntg_id << "\n";
      for (ul i = 0; i < ref_.at(*cntg_id).size(); ++i) {
        fa << cntg_seq_.at(*iso).at(*cntg_id).at(i);
        if (i % 80 == 79) {
          fa << "\n";
        }
      }
      if (ref_.at(*cntg_id).size() % 80 != 0) {
        fa << "\n";
      }
    }
  }
}
