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
#include "sam_l.h"

class Fuga {
 public:
  Fuga(string o_dir, double af, ul dp)
      : kAF(af), kDP(dp), kACGT("ACGT"), totalBase_(0),
        genomeSize_(0), mapped_readNum_(0),
        min_dp_(0), called_posNum_(0),
        AF_fltPos_f_(o_dir + "/allele_freq_filtered_pos"),
        DP_fltPos_f_(o_dir + "/depth_filtered_pos"),
        all_fltPos_f_(o_dir + "/all_filtered_pos"),
        log_f_(o_dir + "/log") {}
  ~Fuga() {}
  void Main(string sam_nom, string ref_nom);

 private:
  const double kAF;
  const ul kDP;
  const string kACGT;
  void ReadRef(string f_nom);
  void ReadSam(string f_nom);
  void SetDP(const SamL& sam_l);
  ul Nuc2Idx(char nuc);
  void Print();
  vector<string> ctgID_;
  unordered_map<string, ul> ctgLen_;
  unordered_map<string, vector<vector<ul>>> dp_;
  ul totalBase_;
  ul genomeSize_;
  ul mapped_readNum_;
  double min_dp_;
  ul called_posNum_;
  ofstream AF_fltPos_f_;
  ofstream DP_fltPos_f_;
  ofstream all_fltPos_f_;
  ofstream log_f_;
};

#endif  // FUGA_H_
