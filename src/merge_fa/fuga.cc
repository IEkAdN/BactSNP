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

void Fuga::Main() {
  SetCntgId();
  ReadFa(kAssFaNom_, &AssFa_);
  ReadFa(kMapFaNom_, &MapFa_);
  PrintConsensusFa();
}

void Fuga::SetCntgId() {
  ifstream AssFaF(kAssFaNom_);
  string AssFaL;
  while (getline(AssFaF, AssFaL)) {
    if (! AssFaL.empty() && AssFaL.at(0) == '>') {
      CntgId_.emplace_back(AssFaL.substr(1, AssFaL.find(" ") - 1));
    }
  }
}

void Fuga::ReadFa(const string& FNom, unordered_map<string, string>* Fa) {
  ifstream F(FNom);
  string L;
  string CntgId;
  while (getline(F, L)) {
    if (! L.empty()) {
      if (L.at(0) == '>') {
        CntgId = L.substr(1, L.find(" ") - 1);
        (*Fa)[CntgId] = "";
      } else {
        Fa->at(CntgId) += L;
      }
    }
  }
}

void Fuga::PrintConsensusFa() {
  for (auto CntgId = CntgId_.begin(); CntgId != CntgId_.end(); ++CntgId) {
    cout << ">" << *CntgId << "\n";
    for (size_t Pos = 0; Pos < AssFa_.at(*CntgId).size(); ++Pos) {
      if (AssFa_.at(*CntgId).at(Pos) ==
          MapFa_.at(*CntgId).at(Pos)) {
        cout << AssFa_.at(*CntgId).at(Pos);
      } else {
        cout << "-";
      }
      if (Pos % 80 == 79) {
        cout << "\n";
      }
    }
    if (AssFa_.at(*CntgId).size() % 80 != 0) {
      cout << "\n";
    }
  }
}
