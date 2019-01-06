#include "hoge.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using std::accumulate;
using std::cerr;
using std::cin;
using std::cout;
using std::distance;
using std::greater;
using std::ifstream;
using std::istream;
using std::make_pair;
using std::map;
using std::max;
using std::min;
using std::multimap;
using std::ofstream;
using std::pair;
using std::stod;
using std::stol;
using std::stoll;
using std::stoul;
using std::stoull;
using std::string;
using std::to_string;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using boost::algorithm::is_any_of;
using boost::algorithm::split;
using boost::algorithm::token_compress_on;

typedef unsigned long long ull;
typedef unsigned long ul;

void Hoge::Main(const string& f_nom) {
  ReadVcf(f_nom);
}

void Hoge::ReadVcf(const string& f_nom) {
  ifstream f(f_nom);
  string l = "";
  vector<string> l_sp;
  while (getline(f, l)) {
    if (l.substr(0, 6) == "#CHROM") {
      split(l_sp, l, is_any_of("\t"));
      cout << "contig\tpos";
      for (ul i = 0; i + 9 < l_sp.size(); ++i) {
        cout << '\t' << l_sp.at(i + 9);
      }
      cout << '\n';
    } else if (l.at(0) != '#') {
      split(l_sp, l, is_any_of("\t"));
      filter_ = l_sp.at(6);
      info_ = l_sp.at(7);
      split(info_sp_, info_, is_any_of(";"));
      for (auto it = info_sp_.begin(); it != info_sp_.end(); ++it) {
        if ((*it).find("SVTYPE") != string::npos) {
          svtype_ = (*it).substr((*it).find('=') + 1, 3);
        }
      }
      if (filter_ == "PASS" && svtype_ == "SNP") {
        Filter(l_sp);
      }
    }
    ResetVar();
  }
}

void Hoge::Filter(const vector<string>& l_sp) {
  string format_(l_sp.at(8));
  if (format_.find("GT") == string::npos ||
      format_.find("COV") == string::npos) {
    return;
  }
  split(format_sp_, format_, is_any_of(":"));
  for (auto it = format_sp_.begin(); it != format_sp_.end(); ++it) {
    if (*it == "GT") {
      gt_idx_ = distance(format_sp_.begin(), it);
    } else if (*it == "COV") {
      cov_idx_ = distance(format_sp_.begin(), it);
    }
  }
  iso_num_ = l_sp.size() - 9;
  cntg_id_ = l_sp.at(0);
  pos_ = stoul(l_sp.at(1)) - 1;
  if (l_sp.at(3).size() != 1 || l_sp.at(4).size() != 1) {
    cerr << "unexpected vcf (1)\n";
    exit(1);
  } else {
    ref_nuc_ = l_sp.at(3).at(0);
    alt_nuc_ = l_sp.at(4).at(0);;
  }
  for (ul i = 0; i < iso_num_; ++i) {
    ReadFormatVal(l_sp, i);
  }
  if (gt_set_.size() != 1) {
    cout << cntg_id_ << '\t' << pos_ + 1;
    for (ul i = 0; i < iso_num_; ++i) {
      cout << '\t' << iso_nuc_.at(i);
    }
    cout << '\n';
  }
}

void Hoge::ReadFormatVal(const vector<string>& l_sp, ul iso_idx) {
  string formatVal(l_sp.at(9 + iso_idx));
  vector<string> formatVal_sp;
  split(formatVal_sp, formatVal, is_any_of(":"));
  ul gt = GetGT(formatVal_sp);
  if (gt != ULONG_MAX) {
    gt_set_.insert(gt);
  }
  SetIsoNuc(gt);
}

ul Hoge::GetGT(const vector<string>& formatVal_sp) {
  string gtStr(formatVal_sp.at(gt_idx_));
  vector<string> gtStr_sp;
  split(gtStr_sp, gtStr, is_any_of("/"));
  if (gtStr_sp.at(0) != gtStr_sp.at(1)) {
    cerr << "This vcf contains hetero GTs\n";
    exit(1);
  }
  if (gtStr_sp.at(0) == ".") {
    return ULONG_MAX;
  } else {
    return stoul(gtStr_sp.at(0));
  }
}


void Hoge::SetIsoNuc(ul gt) {
  if (gt == ULONG_MAX) {
    iso_nuc_.push_back('n');
  } else if (gt == 0) {
    iso_nuc_.push_back(ref_nuc_);
  } else if (gt == 1) {
    iso_nuc_.push_back(alt_nuc_);
  } else {
    cerr << "unexpected vcf (2)\n";
    exit(1);
  }
}

void Hoge::ResetVar() {
  gt_set_.clear();
  iso_nuc_.clear();
}
