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

#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga() {}
  ~Fuga() {}
  void Main(string ref_nom, string region_nom);

 private:
  void ReadFa(string f_nom);
  void ReadRegion(string f_nom);
  void Print();
  unordered_map<string, string> fa_;
  vector<string> fa_header_;
  unordered_map<string, vector<bool> > region_;
};

#endif  // FUGA_H_
