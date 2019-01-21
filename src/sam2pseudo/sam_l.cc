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

# include "sam_l.h"

bool SamL::IsValid() {
  split(l_sp_, l_, "\t");
  flag_ = stoul(l_sp_.at(1));
  if (! (flag_ & 0x4) && ! (flag_ & 0x100)) {
    return true;
  } else {
    return false;
  }
}

void SamL::ReadL() {
  ctgID_ = l_sp_.at(2);
  pos_ = stoul(l_sp_.at(3)) - 1;
  cigar_ = l_sp_.at(5);
  if (cigar_.find_first_not_of("0123456789MIDSH") != string::npos) {
    cerr << "CIGAR contains not-MIDSH characters.\n";
    exit(1);
  }
  seq_ = l_sp_.at(9);
  SplitCigar();
}

void SamL::SplitCigar() {
  split(cigar_op_, cigar_, "0123456789", true);
  vector<string> cigar_num_str;
  split(cigar_num_str, cigar_, "MIDSH", true);
  for (unsigned i = 0; i < cigar_num_str.size(); ++i) {
    cigar_num_.push_back(stoul(cigar_num_str.at(i)));
  }
}
