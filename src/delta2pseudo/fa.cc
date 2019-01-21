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

#include "fa.h"

#include "hoge.h"

void Fa::ReadFa(const string& fa_nom) {
  ifstream fa_f(fa_nom);
  string fa_l = "";
  string id = "";
  while (getline(fa_f, fa_l)) {
    if (fa_l.empty()) {
      continue;
    } else {
      if (fa_l.at(0) == '>') {
        vector<string> fa_l_sp;
        split(fa_l_sp, fa_l, " ");
        id = fa_l_sp.at(0).substr(1);
        id_.push_back(id);
        id_2_seq_.insert(make_pair(id, ""));
      } else {
        transform(fa_l.begin(), fa_l.end(), fa_l.begin(), ::toupper);
        for (ul i = 0; i < fa_l.size(); ++i) {
          if (kACGT.find(fa_l.at(i)) == string::npos) {
            fa_l.replace(i, 1, 1, '-');
          }
        }
        id_2_seq_.at(id) += fa_l;
      }
    }
  }
}
