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
      iso_num_ = l_sp.size() - 9;
      cout << "contig\tpos";
      for (ul i = 0; i + 9 < l_sp.size(); ++i) {
        cout << '\t' << l_sp.at(i + 9);
      }
      cout << '\n';
    } else if (l.at(0) != '#') {
      split(l_sp, l, is_any_of("\t"));
      ReadLine(l_sp);
    }
    ResetVar();
  }
}

void Hoge::ReadLine(const vector<string>& l_sp) {
  string format_(l_sp.at(8));
  split(format_sp_, format_, is_any_of(":"));
  for (auto it = format_sp_.begin(); it != format_sp_.end(); ++it) {
    if (*it == "GT") {
      gt_idx_ = distance(format_sp_.begin(), it);
    }
  }
  cntg_id_ = l_sp.at(0);
  pos_ = stoul(l_sp.at(1)) - 1;
  ref_nuc_ = l_sp.at(3).at(0);
  if (l_sp.at(3).size() == 1) {
    string alt_nuc_str = l_sp.at(4);
    split(alt_nuc_, alt_nuc_str, is_any_of(","));
    for (ul i = 0; i < iso_num_; ++i) {
      ReadFormatVal(l_sp, i);
    }
    unordered_set<char> iso_nuc_set;
    for (ul i = 0; i < iso_num_; ++i) {
      iso_nuc_set.insert(iso_nuc_.at(i));
    }
    iso_nuc_set.erase('n');
    iso_nuc_set.erase('N');
    iso_nuc_set.erase('*');
    if (iso_nuc_set.size() >= 2) {
      cout << cntg_id_ << '\t' << pos_ + 1;
      for (ul i = 0; i < iso_num_; ++i) {
        cout << '\t' << iso_nuc_.at(i);
      }
      cout << '\n';
    }
  }
}

void Hoge::ReadFormatVal(const vector<string>& l_sp, ul iso_idx) {
  string formatVal(l_sp.at(9 + iso_idx));
  vector<string> formatVal_sp;
  split(formatVal_sp, formatVal, is_any_of(":"));
  SetIsoNuc(GetGT(formatVal_sp));
}

ul Hoge::GetGT(const vector<string>& formatVal_sp) {
  string gtStr(formatVal_sp.at(gt_idx_));
  if (gtStr == ".") {
    return ULONG_MAX;
  } else {
    return stoul(gtStr);
  }
}


void Hoge::SetIsoNuc(ul gt) {
  if (gt == 0) {
    iso_nuc_.push_back(ref_nuc_);
  } else if (gt == ULONG_MAX) {
    iso_nuc_.push_back('n');
  } else {
    if (alt_nuc_.at(gt - 1).size() != 1) {
      iso_nuc_.push_back('n');
    } else {
      iso_nuc_.push_back(alt_nuc_.at(gt - 1).at(0));
    }
  }
}

void Hoge::ResetVar() {
  iso_nuc_.clear();
}
