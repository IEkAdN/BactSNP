#include "two_matrix.h"
#include <fstream>
#include <iostream>
#include <regex>
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using std::cout;
using std::ifstream;
using std::map;
using std::regex;
using std::regex_match;
using std::unordered_map;
using std::unordered_set;
using std::string;
using std::vector;
using boost::algorithm::is_any_of;
using boost::algorithm::split;

typedef unsigned long long ull;
typedef unsigned long ul;

vector<string> TwoMatrix::iso_id_;
map<string, string> TwoMatrix::contig_2_ref_seq;
map<string, map<string, string>> TwoMatrix::iso_id_2_contig_2_seq;

void TwoMatrix::ReadRef(const string& ref_na) {
  ifstream ref_f(ref_na);
  string ref_l = "";
  string contig = "";
  while (getline(ref_f, ref_l)) {
    if (!ref_l.empty()) {
      if (ref_l.at(0) == '>') {
        vector<string> ref_l_sp;
        split(ref_l_sp, ref_l, is_any_of(" "));
        contig = ref_l_sp.at(0).substr(1);
        contig_2_ref_seq[contig] = "";
      } else {
        contig_2_ref_seq[contig] += ref_l;
      }
    }
  }
}

void TwoMatrix::ReadVcfLst(const string& vcf_lst_na) {
  ifstream vcf_lst_f(vcf_lst_na);
  string vcf_lst_l = "";
  while (getline(vcf_lst_f, vcf_lst_l)) {
    vector<string> vcf_lst_l_sp;
    split(vcf_lst_l_sp, vcf_lst_l, is_any_of("\t"));
    const string& iso_id = vcf_lst_l_sp.at(0);
    const string& vcf_na = vcf_lst_l_sp.at(1);
    iso_id_.push_back(iso_id);
    ReadVcf(iso_id, vcf_na);
  }
}

void TwoMatrix::ReadVcf(const string& iso_id, const string& vcf_na) {
  unordered_map<string, vector<bool>*> is_called;
  for (auto it = contig_2_ref_seq.begin(); it != contig_2_ref_seq.end();
       ++it) {
    is_called.insert(make_pair((*it).first,
                               new vector<bool>((*it).second.size(), false)));
  }
  iso_id_2_contig_2_seq[iso_id] = contig_2_ref_seq;
  ifstream vcf_f(vcf_na);
  string vcf_l = "";
  while (getline(vcf_f, vcf_l)) {
    if (vcf_l.at(0) != '#' && vcf_l.find("INDEL") == string::npos) {
      vector<string> vcf_l_sp;
      split(vcf_l_sp, vcf_l, is_any_of("\t"));

      const string& contig = vcf_l_sp.at(0);
      const ul pos = stoul(vcf_l_sp.at(1)) - 1;
      (*(is_called.at(contig))).at(pos) = true;

      const string& fmt = vcf_l_sp.at(8);
      const string& fmt_val = vcf_l_sp.at(9);
      vector<string> fmt_sp;
      vector<string> fmt_val_sp;
      split(fmt_sp, fmt, is_any_of(":"));
      split(fmt_val_sp, fmt_val, is_any_of(":"));
      auto gt_itr = find(fmt_sp.begin(), fmt_sp.end(), "GT");
      ul gt_idx = static_cast<ul>(distance(fmt_sp.begin(), gt_itr));
      const string& gt = fmt_val_sp.at(gt_idx);

      const string& info = vcf_l_sp.at(7);
      vector<string> info_sp;
      split(info_sp, info, is_any_of(";"));
      double ref_read_num = 0;
      double alt_read_num = 0;
      for (auto it = info_sp.begin(); it != info_sp.end(); ++it) {
        if ((*it).find("DP4") != string::npos) {
          const string& dp4 = (*it).substr(4);
          vector<string> dp4_sp;
          split(dp4_sp, dp4, is_any_of(","));
          ref_read_num = stod(dp4_sp.at(0)) + stod(dp4_sp.at(1));
          alt_read_num = stod(dp4_sp.at(2)) + stod(dp4_sp.at(3));
        }
      }

      double af = 0;
      string gt_nuc("");
      if (gt == ".") {
        af = 0;
      } else if (gt == "0") {
        gt_nuc = vcf_l_sp.at(3);
        af = ref_read_num / (ref_read_num + alt_read_num);
      } else {
        const string& alt = vcf_l_sp.at(4);
        vector<string> alt_sp;
        split(alt_sp, alt, is_any_of(","));
        gt_nuc = alt_sp.at(stoul(gt) - 1);
        af = alt_read_num / (ref_read_num + alt_read_num);
      }

      const ul qual = stoul(vcf_l_sp.at(5));
      if (qual >= 30 && af >= 0.75) {
        iso_id_2_contig_2_seq[iso_id][contig].replace(pos, 1, gt_nuc);
      } else {
        iso_id_2_contig_2_seq[iso_id][contig].replace(pos, 1, "n");
      }
    }
  }
  for (auto it = contig_2_ref_seq.begin(); it != contig_2_ref_seq.end();
       ++it) {
    for (ul pos = 0; pos < (*it).second.size(); ++pos) {
      if ((*is_called.at((*it).first)).at(pos) == false) {
        iso_id_2_contig_2_seq[iso_id][(*it).first].replace(pos, 1, "n");
      }
    }
  }
  for (auto it = is_called.begin(); it != is_called.end(); ++it) {
    delete (*it).second;
  }
}

void TwoMatrix::PrintMatrix() {
  cout << "contig\tpos";
  for (auto i = iso_id_.begin(); i != iso_id_.end(); ++i) {
    cout << '\t' << (*i);
  }
  cout << '\n';
  for (auto i = contig_2_ref_seq.begin(); i != contig_2_ref_seq.end(); ++i) {
    const string& contig = (*i).first;
    const string& ref_seq = (*i).second;
    for (int64_t pos = 0; pos < static_cast<int64_t>((*i).second.size());
         ++pos) {
      char ref_nuc = ref_seq.at(pos);
      unordered_set<char> alt_nuc_set;
      for (auto j = iso_id_2_contig_2_seq.begin();
          j != iso_id_2_contig_2_seq.end(); ++j) {
        const string& alt_seq = (*j).second[contig];
        char alt_nuc = alt_seq.at(pos);
        if (alt_nuc != 'n' && alt_nuc != 'N') {
          alt_nuc_set.insert(alt_nuc);
        }
      }
      if (alt_nuc_set.size() >= 2) {
        PrintMatrixLine(contig, pos, ref_nuc);
      }
    }
  }
}

void TwoMatrix::PrintMatrixLine(const string& contig, int64_t pos,
                                char ref_nuc) {
    cout << contig << '\t' << pos + 1;
    for (auto i = iso_id_.begin(); i != iso_id_.end(); ++i) {
      const string& alt_seq = iso_id_2_contig_2_seq.at((*i)).at(contig);
      char alt_nuc = alt_seq.at(pos);
      cout << '\t' << alt_nuc;
    }
    cout << '\n';
}

void TwoMatrix::main(const string& vcf_lst_na, const string& ref_na) {
  ReadRef(ref_na);
  ReadVcfLst(vcf_lst_na);
  PrintMatrix();
}
