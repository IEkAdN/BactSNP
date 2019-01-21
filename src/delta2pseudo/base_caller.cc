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

#include "base_caller.h"

void BaseCaller::Main(const string& ref_nom, const string& qry_nom,
                      const string& delta_nom) {
  ReadFa(ref_nom, qry_nom);
  ReadDelta(delta_nom);
  CallBase();
}

void BaseCaller::ReadFa(const string& ref_nom, const string& qry_nom) {
  ref_fa_.ReadFa(ref_nom);
  qry_fa_.ReadFa(qry_nom);
  ref_id_ = ref_fa_.id();
  const unordered_map<string, string>& ref_id_2_seq = ref_fa_.id_2_seq();
  const unordered_map<string, string>& qry_id_2_seq = qry_fa_.id_2_seq();
  for (auto it = ref_id_2_seq.begin(); it != ref_id_2_seq.end(); ++it) {
    ref_id_2_pos_info_[(*it).first].resize((*it).second.size());
  }
  for (auto it = qry_id_2_seq.begin(); it != qry_id_2_seq.end(); ++it) {
    qry_id_2_pos_info_[(*it).first].resize((*it).second.size());
  }
}

void BaseCaller::ReadDelta(const string& delta_nom) {
  ifstream delta_f(delta_nom);
  string delta_l("");
  int aln_idx = 0;
  string ref_id;
  string qry_id;
  getline(delta_f, delta_l);
  getline(delta_f, delta_l);
  while (getline(delta_f, delta_l)) {
    vector<string> delta_l_sp;
    split(delta_l_sp, delta_l, " ");
    if (delta_l.at(0) == '>') {
      ref_id = delta_l_sp.at(0).substr(1);
      qry_id = delta_l_sp.at(1);
    } else if (delta_l_sp.size() == 7) {
      const ll ref_beg = stoll(delta_l_sp.at(0)) - 1;
      const ll ref_end = stoll(delta_l_sp.at(1)) - 1;
      const ll qry_beg = stoll(delta_l_sp.at(2)) - 1;
      const ll qry_end = stoll(delta_l_sp.at(3)) - 1;
      add_2_aln(ref_beg, qry_beg, qry_end, ref_id, qry_id);
      string aln_ref_seq = 
        ref_fa_.seq(ref_id).substr(ref_beg, ref_end - ref_beg + 1);
      aln_.at(aln_idx).set_ref_seq(aln_ref_seq);
      for (ll i = ref_beg; i != ref_end + 1; ++i) {
        ref_id_2_pos_info_.at(ref_id).at(i).SetInfo(aln_idx);
      }
      if (qry_beg < qry_end) {
        for (ll i = qry_beg; i != qry_end + 1; ++i) {
          qry_id_2_pos_info_.at(qry_id).at(i).SetInfo('+');
        }
      } else {
        for (ll i = qry_end; i != qry_beg + 1; ++i) {
          qry_id_2_pos_info_[qry_id].at(i).SetInfo('-');
        }
      }
      ++aln_idx;
    } else {
      aln_.at(aln_idx - 1).add_2_digit(stoi(delta_l));
    }
  }
  for (auto i = aln_.begin(); i != aln_.end(); ++i) {
    (*i).ReadDigit();
  }
}

void BaseCaller::CallBase() const {
  for (auto i = ref_id_.begin(); i != ref_id_.end(); ++i) {
    const string& ref_id = *i;
    cout << '>' << ref_id << '\n';
    const vector<RefPosInfo>& ref_pos_info_v = ref_id_2_pos_info_.at(ref_id);
    for (long long j = 0; 
         j != static_cast<long long>(ref_pos_info_v.size()); ++j) {
      const vector<int>& aln_idx = ref_pos_info_v.at(j).aln_idx();
      unordered_set<char> nuc;
      for (auto k = aln_idx.begin(); k != aln_idx.end(); ++k) {
        Aln aln(aln_.at(*k));
        const string& qry_id = aln.qry_id();
        const long long qry_pos = aln.RefPos2QryPos(j);
        if (! aln.is_near_indel(j) && qry_pos != -1) {
          if (qry_id_2_pos_info_.at(qry_id).at(qry_pos).strand() == '+') {
            nuc.insert(qry_fa_.nuc(qry_id, qry_pos));
          } else {
            nuc.insert(GetComplement(qry_fa_.nuc(qry_id, qry_pos)));
          }
        } else {
          nuc.insert('-');
        }
      }
      if (nuc.find('-') != nuc.end() || nuc.size() != 1) {
        cout << '-';
      } else {
        cout << *nuc.begin();
      }
      if (j % 80 == 79) {
        cout << '\n';
      }
    }
    if (ref_pos_info_v.size() % 80 != 0) {
      cout << '\n';
    }
  }
}

char BaseCaller::GetComplement(const char nuc) const {
  if (nuc == 'A') {
    return 'T';
  } else if (nuc == 'C') {
    return 'G';
  } else if (nuc == 'G') {
    return 'C';
  } else if (nuc == 'T') {
    return 'A';
  } else {
    return '-';
  }
}

void Aln::ReadDigit() {
  ref_gap_.resize(ref_seq_.size(), 0);
  ref_near_indel_.resize(ref_seq_.size(), false);
  int ref_del_num = 0;
  for (auto i = digit_.begin(); i != digit_.end(); ++i) {
    if ((*i) < 0) {
      ++ref_del_num;
    }
  }
  qry_gap_.resize(ref_seq_.size() + ref_del_num, 0);

  ll done_ref_len = 0;
  ll digit_sum = 0;
  for (auto i = digit_.begin(); i != digit_.end(); ++i) {
    digit_sum += abs(*i);
    if ((*i) > 0) {
      done_ref_len += (*i);
      qry_gap_.at(digit_sum - 1) = INT_MIN;
      for (ll j = digit_sum; j < static_cast<ll>(qry_gap_.size()); ++j) {
        --(qry_gap_.at(j));
      }
      for (ll j = max(static_cast<long long>(0), done_ref_len - indel_);
           j <= min(static_cast<long long>(ref_seq_.size() - 1),
                    done_ref_len + indel_);
           ++j) {
        ref_near_indel_.at(j) = true;
      }
    } else if ((*i) < 0) {
      done_ref_len += -(*i) - 1;
      for (ll j = done_ref_len; j < static_cast<ll>(ref_gap_.size()); ++j) {
        ++(ref_gap_.at(j));
      }
      for (ll j = max(static_cast<long long>(0), done_ref_len + 1 - indel_);
           j <= min(static_cast<long long>(ref_seq_.size() - 1),
                    done_ref_len + indel_);
           ++j) {
        ref_near_indel_.at(j) = true;
      }
    }
  }
}

ll Aln::RefPos2QryPos(const ll ref_pos) const {
  const ll aln_pos = RefPos2AlnPos(ref_pos);
  return AlnPos2QryPos(aln_pos);
}

ll Aln::RefPos2AlnPos(const ll ref_pos) const {
  ll aln_ref_pos = ref_pos - ref_beg_;
  return aln_ref_pos + ref_gap_.at(aln_ref_pos);
}

ll Aln::AlnPos2QryPos(const ll aln_pos) const {
  if (qry_gap_.at(aln_pos) == INT_MIN) {
    return -1;
  } else { 
    ll aln_qry_pos = aln_pos + qry_gap_.at(aln_pos);
    if (qry_beg_ < qry_end_) {
      return qry_beg_ + aln_qry_pos;
    } else {
      return qry_beg_ - aln_qry_pos;
    }
  }
}
