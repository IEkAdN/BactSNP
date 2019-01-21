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

#ifndef SAM_L_H_
#define SAM_L_H_

#include "hoge.h"

class SamL {
 public:
  SamL(const string& l) : l_(l), flag_(0), pos_(0) {}
  ~SamL() {}
  bool IsValid();
  void ReadL();
  ul pos() const { return pos_; }
  const string& ctgID() const { return ctgID_; }
  const string& seq() const { return seq_; }
  const vector<string> cigar_op() const { return cigar_op_; }
  const vector<ul> cigar_num() const { return cigar_num_; }

 private:
  void SplitCigar();
  string l_;
  vector<string> l_sp_;
  ul flag_;
  string ctgID_;
  ul pos_;
  string seq_;
  vector<string> cigar_op_;
  vector<ul> cigar_num_;
  string cigar_;
};

#endif  // SAM_L_H_
