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
    split(l_sp, l, is_any_of("\t"));
    if (l.substr(0, 6) == "#CHROM") {
      cout << "contig\tpos";
      for (ul i = 0; i + 9 < l_sp.size(); ++i) {
        cout << '\t' << l_sp.at(i + 9);
      }
      cout << '\n';
    } else if (l.at(0) != '#') {
      if (! LineIsExpectedFmt(l, l_sp)) {
        exit(1);
      } else {
        PrintInterIsoSNP(l, l_sp);
      }
    }
  }
}

bool Hoge::LineIsExpectedFmt(const string& l, const vector<string>& l_sp) {
  string format(l_sp.at(8));
  if (format.find("GT") == string::npos) {
    cerr << "the line below doesn't contain GT\n"
         << l << '\n';
    return false;
  }
  ul type_idx = l.find("TYPE=");
  ul semicolon_idx = l.find_first_of(";\t", type_idx);
  string type = l.substr(type_idx + 5, semicolon_idx - type_idx - 5);
  vector<string> type_sp;
  split(type_sp, type, is_any_of(","));
  for (auto it = type_sp.begin(); it != type_sp.end(); ++it) {
    if ((*it) != "snp" && (*it) != "mnp") {
      cerr << "the line below contains non-snp or non-mnp TYPE call\n"
           << l << '\n';
      return false;
    }
  }
  string ref_nuc = l_sp.at(3);
  string alt_nuc = l_sp.at(4);
  vector<string> alt_nuc_sp;
  split(alt_nuc_sp, alt_nuc, is_any_of(","));
  for (auto it = alt_nuc_sp.begin(); it != alt_nuc_sp.end(); ++it) {
    if ((*it).size() != ref_nuc.size()) {
      cerr << "the line below contains ALT whose len. differs from REF len.\n"
           << l << '\n';
      return false;
    }
  }
  return true;
}

void Hoge::PrintInterIsoSNP(const string& l, const vector<string>& l_sp) {
  string format(l_sp.at(8));
  vector<string> format_sp;
  split(format_sp, format, is_any_of(":"));
  ul gt_idx(0);
  for (auto it = format_sp.begin(); it != format_sp.end(); ++it) {
    if (*it == "GT") {
      gt_idx = distance(format_sp.begin(), it);
    }
  }
  ul iso_num = l_sp.size() - 9;
  string cntg_id = l_sp.at(0);
  ul pos = stoul(l_sp.at(1)) - 1;
  string ref_nuc = l_sp.at(3);
  string alt_nuc = l_sp.at(4);
  vector<string> alt_nuc_sp;
  split(alt_nuc_sp, alt_nuc, is_any_of(","));
  vector<ul> gtVec;
  for (ul iso_idx = 0; iso_idx < iso_num; ++iso_idx) {
    gtVec.push_back(GetGT(l_sp, gt_idx, iso_idx));
  }
  for (ul pos_idx = 0; pos_idx < ref_nuc.size(); ++pos_idx) {
    vector<char> isoNuc;
    for (ul iso_idx = 0; iso_idx < iso_num; ++iso_idx) {
      ul gt = gtVec.at(iso_idx);
      if (gt == 0) {
        isoNuc.push_back(ref_nuc.at(pos_idx));
      } else if (gt == ULONG_MAX) {
        isoNuc.push_back('n');
      } else {
        isoNuc.push_back(alt_nuc_sp.at(gt - 1).at(pos_idx));
      }
    }
    unordered_set<char> isoNucSet(isoNuc.begin(), isoNuc.end());
    isoNucSet.erase('n');
    if (isoNucSet.size() != 1) {
      cout << cntg_id << '\t' << pos + 1 + pos_idx;
      for (ul i = 0; i < iso_num; ++i) {
        cout << '\t' << isoNuc.at(i);
      }
      cout << '\n';
    }
  }
}

ul Hoge::GetGT(const vector<string>& l_sp, ul gt_idx, ul iso_idx) {
  string formatVal(l_sp.at(9 + iso_idx));
  vector<string> formatVal_sp;
  split(formatVal_sp, formatVal, is_any_of(":"));
  string gt(formatVal_sp.at(gt_idx));
  if (gt == ".") {
    return ULONG_MAX;
  } else {
    return stoul(gt);
  }
}
