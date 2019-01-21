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

#ifndef BASE_CALLER_H
#define BASE_CALLER_H

#include "base_caller.h"

#include <climits>

#include "hoge.h"
#include "fa.h"

class BaseCaller;
class Aln;
class RefPosInfo;
class QryPosInfo;

class Aln {
 public:
  Aln(const long long ref_beg, const long long qry_beg, const long long qry_end,
      const std::string& ref_id, const std::string& qry_id,
      const long long indel)
      : ref_beg_(ref_beg),
        qry_beg_(qry_beg),
        qry_end_(qry_end),
        ref_id_(ref_id),
        qry_id_(qry_id),
        indel_(indel) {}
  ~Aln() {}
  void set_ref_seq(const std::string& aln_ref_seq) {
    ref_seq_ = aln_ref_seq;
  }
  void add_2_digit(const int digit) { digit_.push_back(digit); }
  bool is_near_indel(const long long ref_pos) const {
    return ref_near_indel_.at(ref_pos - ref_beg_);
  }
  void ReadDigit();
  const std::string& qry_id() const { return qry_id_; }
  long long RefPos2QryPos(const long long ref_pos) const;
  void SetAlnSeq(long long i);
 
 private:
  long long ref_beg_;
  long long qry_beg_;
  long long qry_end_;
  std::string ref_id_;
  std::string qry_id_;
  long long indel_;
  std::string ref_seq_;
  std::vector<int> digit_;
  std::vector<int> ref_gap_;
  std::vector<int> qry_gap_;
  std::vector<bool> ref_near_indel_;
  long long RefPos2AlnPos(const long long ref_pos) const;
  long long AlnPos2QryPos(const long long aln_pos) const;
};

class BaseCaller {
 public:
  const long long kIndel;
  BaseCaller(const long long indel) : kIndel(indel) {}
  ~BaseCaller() {}
  void Main(const std::string& ref_nom, const std::string& qry_nom,
            const std::string& delta_nom);

 private:
  void ReadFa(const std::string& ref_nom, const std::string& qry_nom);
  void ReadDelta(const std::string& delta_nom);
  void add_2_aln(const long long ref_beg, const long long qry_beg,
                 const long long qry_end, const std::string& ref_id,
                 const std::string& qry_id) {
    aln_.push_back(Aln(ref_beg, qry_beg, qry_end, ref_id, qry_id, kIndel));
  }
  void CallBase() const;
  char GetComplement(const char nuc) const;
  Fa ref_fa_;
  Fa qry_fa_;
  std::vector<Aln> aln_;
  std::vector<std::string> ref_id_;
  std::map<std::string, std::vector<RefPosInfo>> ref_id_2_pos_info_;
  std::map<std::string, std::vector<QryPosInfo>> qry_id_2_pos_info_;
};

class RefPosInfo {
 public:
  RefPosInfo() : depth_(0) {}
  ~RefPosInfo() {}
  int depth() const { return depth_; }
  vector<int> aln_idx() const { return aln_idx_; }
  void SetInfo(const long long aln_idx) {
    aln_idx_.emplace_back(aln_idx);
    ++depth_;
  }

 private:
  int depth_;
  vector<int> aln_idx_;
};

class QryPosInfo {
 public:
  QryPosInfo()
      : depth_(0) {}
  ~QryPosInfo() {}
  int depth() const { return depth_; }
  char strand() const { return strand_; }
  void SetInfo(char strand) {
    strand_ = strand;
    ++depth_;
  }

 private:
  int depth_;
  char strand_;
};

#endif  // BASE_CALLER_H
