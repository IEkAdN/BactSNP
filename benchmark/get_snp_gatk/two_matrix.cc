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

void TwoMatrix::ReadVcf(const string& vcf_na) {
  ifstream vcf_f(vcf_na);
  string vcf_l = "";
  while (getline(vcf_f, vcf_l)) {
    vector<string> vcf_l_sp;
    split(vcf_l_sp, vcf_l, is_any_of("\t"));
    if (vcf_l.at(0) == '#') {
      if (vcf_l.find("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
          != string::npos) {
        cout << "contig\tpos";
        for (ul i = 9; i < vcf_l_sp.size(); ++i) {
          iso_id_.push_back(vcf_l_sp.at(i));
          cout << '\t' << vcf_l_sp.at(i);
        }
        cout << '\n';
      }
    } else {
      const string& contig = vcf_l_sp.at(0);
      const ul pos = stoul(vcf_l_sp.at(1)) - 1;
      const string& ref = vcf_l_sp.at(3);
      const string& alt = vcf_l_sp.at(4);
      const string& filter = vcf_l_sp.at(6);
      if (filter == "PASS") {
        const string& fmt = vcf_l_sp.at(8);
        vector<string> fmt_sp;
        split(fmt_sp, fmt, is_any_of(":"));
        auto gt_itr = find(fmt_sp.begin(), fmt_sp.end(), "GT");
        ul gt_idx = static_cast<ul>(distance(fmt_sp.begin(), gt_itr));
        vector<string> gt(iso_id_.size());
        for (ul iso_idx = 0; iso_idx < iso_id_.size(); ++iso_idx) {
          vector<string> fmt_val_sp;
          split(fmt_val_sp, vcf_l_sp.at(9 + iso_idx), is_any_of(":"));
          gt.at(iso_idx) = fmt_val_sp.at(gt_idx);
        }
        vector<string> alt_sp;
        split(alt_sp, alt, is_any_of(","));
        unordered_set<string> gt_set;
        for (auto i = gt.begin(); i != gt.end(); ++i) {
          if (*i == "0" || (*i != "." && alt_sp.at(stoul(*i) - 1) != "*")) {
            gt_set.insert(*i);
          }
        }
        if (gt_set.size() >= 2) {
          cout << contig << '\t' << pos + 1;
          for (ul iso_idx = 0; iso_idx < iso_id_.size(); ++iso_idx) {
            if (gt.at(iso_idx) == "0") {
              cout << '\t' << ref;
            } else if (gt.at(iso_idx) == "." || alt_sp.at(stoul(gt.at(iso_idx)) - 1) == "*") {
              cout << '\t' << 'n';
            } else {
              cout << '\t' << alt_sp.at(stoul(gt.at(iso_idx)) - 1);
            }
          }
          cout << '\n';
        }
      }
    }
  }
}

void TwoMatrix::main(const string& vcf_na) {
  ReadVcf(vcf_na);
}
