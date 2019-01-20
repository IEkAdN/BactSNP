#include "fuga.h"

void Fuga::Main() {
  ReadRoot();
  SetAln();
  SetOrgSnp();
  ReadDelta();
  MoveSnp();
}

void Fuga::ReadRoot() {
  ifstream f(kRootGenomeFNom_);
  string l = "";
  vector<string> l_sp;
  string root_cntg("");
  while (getline(f, l)) {
    if (! l.empty()) {
      if (l.at(0) == '>') {
        split(l_sp, l, " ");
        root_cntg = l_sp.at(0).substr(1);
        root_cntg_.push_back(root_cntg);
        root_.insert(make_pair(root_cntg, ""));
      } else {
        root_.at(root_cntg) += l;
        root_merged_seq_ += l;
        root_size_ += l.size();
      }
    }
  }
}

void Fuga::SetAln() {
  ifstream F(kRootAlnFNom_);
  string L;
  getline(F, L);
  if (L != ">Node1 The root") {
    cerr << "ERROR: Unexpected fasta\n";
    exit(1);
  }
  getline(F, RootSeq_);
  u32 GapNum(count(L.begin(), L.end(), '-'));
  AlnLen_ = RootSeq_.size();
  AlnPos2Pos_.resize(AlnLen_, 0);
  Pos2AlnPos_.resize(AlnLen_ - GapNum, 0);
  u32 GapCntr(0);
  for (u32 i = 0; i < AlnLen_; ++i) {
    if (RootSeq_.at(i) == '-') {
      ++GapCntr;
      AlnPos2Pos_.at(i) = UINT32_MAX;
    } else {
      AlnPos2Pos_.at(i) = i - GapCntr;
      Pos2AlnPos_.at(i - GapCntr) = i;
    }
  }
}

void Fuga::SetOrgSnp() {
  ifstream F(kIsoAlnFNom_);
  string L;
  while (getline(F, L)) {
    if (L.at(0) == '>') {
      if (L != ">A10") {
        Iso_.emplace_back(L.substr(1, 1) + "000" + L.substr(2, 1));
      } else {
        Iso_.emplace_back("A0010");
      }
    } else {
      IsoSeq_.emplace_back(L);
    }
  }
  IsoNum_ = Iso_.size();
  for (u32 i = 0; i < AlnLen_; ++i) {
    unordered_set<char> Nuc;
    for (auto j = IsoSeq_.begin(); j != IsoSeq_.end(); ++j) {
      char _Nuc(j->at(i));
      if (_Nuc == 'A' || _Nuc == 'C' || _Nuc == 'G' || _Nuc == 'T') {
        Nuc.insert(_Nuc);
      }
    }
    if (Nuc.size() == 2) {
      OrgSnpPos_.emplace_back(i);
    } else if (Nuc.size() >= 3) {
      cerr << "This position has 3 or more alleles:\n" << i << '\n';
      exit(1);
    }
  }
  SnpNum_ = OrgSnpPos_.size();
  RootOrAlt_.resize(SnpNum_);
  for (auto i = RootOrAlt_.begin(); i != RootOrAlt_.end(); ++i) {
    i->resize(IsoNum_);
  }
  for (u32 i = 0; i < SnpNum_; ++i) {
    u32 OrgSnpPos(OrgSnpPos_.at(i));
    char RootNuc(RootSeq_.at(OrgSnpPos));
    for (u32 j = 0; j < IsoNum_; ++j) {
      char Nuc(IsoSeq_.at(j).at(OrgSnpPos));
      if (Nuc == RootNuc) {
        RootOrAlt_.at(i).at(j) = true;
      } else {
        RootOrAlt_.at(i).at(j) = false;
      }
    }
  }
}

void Fuga::ReadDelta() {
  ifstream f(kDeltaFNom_);
  string l = "";
  vector<string> l_sp;
  int aln_idx = 0;
  string ref_id;
  string qry_id;
  getline(f, l);
  getline(f, l);
  while (getline(f, l)) {
    split(l_sp, l, " ");
    if (l.at(0) == '>') {
      ref_id = l_sp.at(0).substr(1);
      qry_id = l_sp.at(1);
      const int64_t ref_len = stoll(l_sp.at(2));
      const int64_t qry_len = stoll(l_sp.at(3));
      if (ref_id_2_pos_info_.find(ref_id) == ref_id_2_pos_info_.end()) {
        ref_id_2_pos_info_[ref_id].resize(ref_len);
      }
      if (qry_id_2_pos_info_.find(qry_id) == qry_id_2_pos_info_.end()) {
        qry_id_2_pos_info_[qry_id].resize(qry_len);
      }
    } else if (l_sp.size() == 7) {
      const int64_t ref_beg = stoll(l_sp.at(0)) - 1;
      const int64_t ref_end = stoll(l_sp.at(1)) - 1;
      const int64_t qry_beg = stoll(l_sp.at(2)) - 1;
      const int64_t qry_end = stoll(l_sp.at(3)) - 1;
      char qry_strand = '0';
      if (qry_beg < qry_end) {
        qry_strand = '+';
      } else {
        qry_strand = '-';
      }
      add_2_aln(ref_beg, qry_beg, ref_end, qry_end, qry_id, qry_strand);
      for (int i = ref_beg; i != ref_end + 1; ++i) {
        ref_id_2_pos_info_.at(ref_id).at(i).SetInfo(aln_idx);
      }
      for (int64_t i = min(qry_beg, qry_end);
           i != max(qry_beg, qry_end) + 1; ++i) {
        qry_id_2_pos_info_[qry_id].at(i).SetInfo();
      }
      ++aln_idx;
    } else {
      aln_.at(aln_.size() - 1).SetSeq(stoll(l));
    }
  }
}

void Fuga::MoveSnp() {
  srand(kRandSeed_);
  cout << "contig\tpos";
  for (auto i = Iso_.begin(); i != Iso_.end(); ++i) {
    cout << "\t" << *i;
  }
  cout << "\n";
  for (u32 i = 0; i < SnpNum_; ++i) {
    pair<string, u32> Tmp = GetUniqPos();
    const string& RootCntg(Tmp.first);
    const u32 RootPos(Tmp.second);
    const u32 AlnIdx(ref_id_2_pos_info_.at(RootCntg).at(RootPos).aln_idx());
    const string& RefCntg(aln_.at(AlnIdx).qry_id());
    const u32 RefPos(aln_.at(AlnIdx).RefPos2QryPos(RootPos));
    char QryStrand(aln_.at(AlnIdx).qry_strand());
    u32 AlnPos(Pos2AlnPos(RootPos));
    char RootNuc(RootSeq_.at(AlnPos));
    char AltRootNuc(GetAltAllele(RootNuc));
    cout << RefCntg << "\t" << RefPos + 1;
    for (u32 j = 0; j < IsoNum_; ++j) {
      IsoSeq_.at(j).replace(OrgSnpPos_.at(i), 1, 1, RootSeq_.at(OrgSnpPos_.at(i)));
      if (RootOrAlt_.at(i).at(j)) {
        if (QryStrand == '+') {
          cout << "\t" <<  RootNuc;
        } else {
          cout << "\t" <<  GetRevComp(RootNuc);
        }
      } else {
        IsoSeq_.at(j).replace(AlnPos, 1, 1, AltRootNuc);
        if (QryStrand == '+') {
          cout << "\t" << AltRootNuc;
        } else {
          cout << "\t" << GetRevComp(AltRootNuc);
        }
      }
    }
    cout << "\n";
  }
  for (u32 i = 0; i < IsoNum_; ++i) {
    ofstream Fa(kOutDir_ + "/" + Iso_.at(i) + ".fa");
    Fa << '>' << root_cntg_.at(0) << "\n";
    for (auto j = IsoSeq_.at(i).begin(); j != IsoSeq_.at(i).end(); ++j) {
      if (*j != '-') {
        Fa << *j;
      }
    }
    Fa << '\n';
  }
}

pair<string, u32> Fuga::GetUniqPos() {
  string cntg("");
  int64_t pos(0);
  while (1) {
    u32 merged_seq_pos = rand() % root_size_;
    u32 cntg_len_sum(0);
    for (auto it = root_cntg_.begin(); it != root_cntg_.end(); ++it) {
      cntg_len_sum += root_.at(*it).size();
      if (cntg_len_sum > merged_seq_pos) {
        cntg = *it;
        pos = merged_seq_pos - (cntg_len_sum - root_.at(*it).size());
        break;
      }
    }
    const RefPosInfo& ref_pos_info = ref_id_2_pos_info_.at(cntg).at(pos);
    if (JdgUniq(cntg, pos, 'r')) {
      const int aln_idx = ref_pos_info.aln_idx();
      const string& qry_id = aln_.at(aln_idx).qry_id();
      const int64_t qry_pos = aln_.at(aln_idx).RefPos2QryPos(pos);
      if (qry_pos != -1) {
        if (JdgUniq(qry_id, qry_pos, 'q')) {
          return make_pair(cntg, pos);
        }
      }
    }
  }
}

bool Fuga::JdgUniq(const string& contig_id, const int64_t pos, const char r_or_q) const {
  int64_t first_ununiq_pos_after_pos = pos;
  int64_t last_ununiq_pos_b4_pos = pos;
  if (r_or_q == 'r') {
    const vector<RefPosInfo>& ref_pos_info = ref_id_2_pos_info_.at(contig_id);
    while (ref_pos_info.at(first_ununiq_pos_after_pos).depth() == 1) {
      if (first_ununiq_pos_after_pos == int64_t(ref_pos_info.size()) - 1) {
        break;
      }
      ++first_ununiq_pos_after_pos;
    }
    while (ref_pos_info.at(last_ununiq_pos_b4_pos).depth() == 1) {
      if (last_ununiq_pos_b4_pos == 0) {
        break;
      }
      --last_ununiq_pos_b4_pos;
    }
    if (first_ununiq_pos_after_pos - last_ununiq_pos_b4_pos - 1 >=
        kUniqAlnLen_) {
      return true;
    } else {
      return false;
    }
  } else {
    const vector<QryPosInfo>& qry_pos_info = qry_id_2_pos_info_.at(contig_id);
    while (qry_pos_info.at(first_ununiq_pos_after_pos).depth() == 1) {
      if (first_ununiq_pos_after_pos == int64_t(qry_pos_info.size()) - 1) {
        break;
      }
      ++first_ununiq_pos_after_pos;
    }
    while (qry_pos_info.at(last_ununiq_pos_b4_pos).depth() == 1) {
      if (last_ununiq_pos_b4_pos == 0) {
        break;
      }
      --last_ununiq_pos_b4_pos;
    }
    if (first_ununiq_pos_after_pos - last_ununiq_pos_b4_pos - 1 >=
        kUniqAlnLen_) {
      return true;
    } else {
      return false;
    }
  }
}

char Fuga::GetAltAllele(char root_allele) {
  const string acgt("ACGT");
  while (1) {
    char var_allele = acgt.at(rand() % 4);
    if (var_allele != root_allele) {
      return var_allele;
    }
  }
}

char Fuga::GetRevComp(char nuc) const {
  if (nuc == 'A') {
    return 'T';
  } else if (nuc == 'C') {
    return 'G';
  } else if (nuc == 'G') {
    return 'C';
  } else if (nuc == 'T') {
    return 'A';
  } else {
    cerr << "unexpected base in the input snp matrix\n";
    exit(1);
  }
}

void Aln::SetSeq(int64_t i) {
  if (i < 0) {
    ref_aln_seq_ += string(abs(i) - 1, 'n') + '.';
    qry_aln_seq_ += string(abs(i), 'n');
  } else if (i > 0) {
    ref_aln_seq_ += string(i, 'n');
    qry_aln_seq_ += string(i - 1, 'n') + '.';
  } else {
    ref_aln_seq_ += string(
        ref_len_ - count(ref_aln_seq_.begin(), ref_aln_seq_.end(), 'n'), 'n');
    qry_aln_seq_ += string(
        qry_len_ - count(qry_aln_seq_.begin(), qry_aln_seq_.end(), 'n'), 'n');
  }
}

int64_t Aln::RefPos2QryPos(const int64_t ref_pos) const {
  const int64_t aln_pos = RefPos2AlnPos(ref_pos);
  return AlnPos2QryPos(aln_pos);
}

int64_t Aln::RefPos2AlnPos(int64_t ref_pos) const {
  auto i = ref_aln_seq_.begin();
  advance(i, ref_pos - ref_beg_);
  while (count(ref_aln_seq_.begin(), i + 1, 'n') < ref_pos - ref_beg_ + 1) {
    ++i;
  }
  return distance(ref_aln_seq_.begin(), i);
}

int64_t Aln::AlnPos2QryPos(int64_t aln_pos) const {
  if (qry_aln_seq_.at(aln_pos) == '.') {
    return -1;
  } else if (qry_strand_ == '+') {
    return
        qry_beg_ +
        count(qry_aln_seq_.begin(), qry_aln_seq_.begin() + aln_pos + 1, 'n') -
        1;
  } else {
    return
        qry_beg_ -
        count(qry_aln_seq_.begin(), qry_aln_seq_.begin() + aln_pos + 1, 'n') +
        1;
  }
}
