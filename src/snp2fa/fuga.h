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
  void Main(string snp_nom, string ref_nom, string out_dir_nom);

 private:
  void ReadRef(string f_nom);
  void ReadSnp(string f_nom);
  bool IsAmbiguous(const vector<string>& l_sp);
  void Print(string dir_nom);
  unordered_map<string, string> ref_;
  vector<string> cntg_id_;
  vector<string> iso_;
  // isolate name -> contig name -> sequence
  unordered_map<string, unordered_map<string, string> > cntg_seq_;
};

#endif  // FUGA_H_
