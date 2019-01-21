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
#include "sam_l.h"

void Fuga::Main(string sam_nom, string ref_nom) {
  ReadRef(ref_nom);
  ReadSam(sam_nom);
  Print();
}

void Fuga::ReadRef(string f_nom) {
  ifstream f(f_nom);
  string l("");
  while (getline(f, l)) {
    if (! l.empty()) {
      if (l.at(0) == '>') {
        ctgID_.push_back(l.substr(1, l.find(" ") - 1));
        ctgLen_[ctgID_.back()] = 0;
      } else {
        ctgLen_.at(ctgID_.back()) += l.size();
        genomeSize_ += l.size();
      }
    }
  }
  for (auto _ctgID = ctgID_.begin(); _ctgID < ctgID_.end(); ++_ctgID) {
    dp_[*_ctgID].resize(ctgLen_.at(*_ctgID), {0, 0, 0, 0, 0, 0});
  }
}

void Fuga::ReadSam(string f_nom) {
  ifstream f(f_nom);
  string l = "";
  std::vector<std::string> l_sp;
  while (getline(f, l)) {
    if (l.at(0) != '@') {
      SamL sam_l(l);
      if (sam_l.IsValid()) {
        ++mapped_readNum_;
        sam_l.ReadL();
        SetDP(sam_l);
      }
    }
  }
}

void Fuga::SetDP(const SamL& sam_l) {
  long cur_ref_pos(long(sam_l.pos()) - 1);
  long cur_read_pos(-1);
  for (ul cigar_idx = 0; cigar_idx < sam_l.cigar_op().size(); ++cigar_idx) {
    string _cigar_op(sam_l.cigar_op().at(cigar_idx));
    ul _cigar_num(sam_l.cigar_num().at(cigar_idx));
    if (_cigar_op == "M") {
      totalBase_ += _cigar_num;
      for (ul i = 0; i < _cigar_num; ++i) {
        ++cur_ref_pos;
        ++cur_read_pos;
        ++(dp_.at(sam_l.ctgID()).at(ul(cur_ref_pos)).at(Nuc2Idx(sam_l.seq().at(ul(cur_read_pos)))));
      }
    } else if (_cigar_op == "D") {
      for (ul i = 0; i < _cigar_num; ++i) {
        ++cur_ref_pos;
        ++(dp_.at(sam_l.ctgID()).at(ul(cur_ref_pos)).at(5));
      }
    } else if (_cigar_op == "S") {
      if (cigar_idx == 0) {
        cur_read_pos += _cigar_num;
      }
    } else if (_cigar_op == "I") {
      cur_read_pos += _cigar_num;
    }
  }
}

ul Fuga::Nuc2Idx(char nuc) {
  switch(nuc) {
    case 'A':
      return 0;
      break;
    case 'C':
      return 1;
      break;
    case 'G':
      return 2;
      break;
    case 'T':
      return 3;
      break;
    case 'N':
      return 4;
      break;
    default:
      log_f_ << "can't convert '" << nuc << "'\n";
      exit(1);
  }
}

void Fuga::Print() {
  min_dp_ = kDP;
  for (auto _ctgID = ctgID_.begin(); _ctgID != ctgID_.end(); ++_ctgID) {
    cout << '>' << *_ctgID << '\n';
    for (ul pos = 0; pos < ctgLen_.at(*_ctgID); ++pos) {
      bool AF_flt(false);
      bool DP_flt(false);
      auto dp_max_it = max_element(dp_.at(*_ctgID).at(pos).begin(),
                                   dp_.at(*_ctgID).at(pos).end() - 2);
      ul dp_sum = accumulate(dp_.at(*_ctgID).at(pos).begin(),
                             dp_.at(*_ctgID).at(pos).end(), 0);
      if (*dp_max_it < dp_sum * kAF) {
        AF_fltPos_f_ << *_ctgID << '\t' << pos + 1 << '\n';
        AF_flt = true;
      }
      if (*dp_max_it < min_dp_) {
        DP_fltPos_f_ << *_ctgID << '\t' << pos + 1 << '\n';
        DP_flt = true;
      }
      if (AF_flt || DP_flt) {
        all_fltPos_f_ << *_ctgID << '\t' << pos + 1 << '\n';
        cout << '-';
      } else {
        ++called_posNum_;
        cout << string("ACGT").at(
            distance(dp_.at(*_ctgID).at(pos).begin(), dp_max_it));
      }
      if (pos % 80 == 79) {
        cout << '\n';
      }
    }
    if (ctgLen_.at(*_ctgID) % 80 != 0) {
      cout << '\n';
    }
  }
  log_f_ << "# mapped read:	" << mapped_readNum_ << '\n';
  log_f_ << "average coverage depth:	" << double(totalBase_) / genomeSize_ << '\n';
  log_f_ << "minimum depth:	" << min_dp_ << '\n';
  log_f_ << "called sites / all sites (%):	"
       << double(called_posNum_) / genomeSize_ * 100 << '\n';
}
