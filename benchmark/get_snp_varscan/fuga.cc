#include "fuga.h"

void Fuga::Main(string result_nom, string isoLst_nom) {
  ReadIsoLst(isoLst_nom);
  ReadResult(result_nom);
}

void Fuga::ReadIsoLst(string f_nom) {
  ifstream f(f_nom);
  string l("");
  vector<string> l_sp;
  cout << "contig\tpos";
  while (getline(f, l)) {
    cout << '\t' << l;
  }
  cout << '\n';
}

void Fuga::ReadResult(string f_nom) {
  ifstream f(f_nom);
  string l("");
  vector<string> l_sp;
  getline(f, l);
  while (getline(f, l)) {
    split(l_sp, l, is_any_of("\t"));
    string ctg_nom(l_sp.at(0));
    ul pos(stoul(l_sp.at(1)) - 1);
    const string& strandFlt(l_sp.at(5).substr(0, 4));
    if (strandFlt != "Pass") {
      continue;
    }
    const string& splInfo(l_sp.at(10));
    vector<string> splInfo_sp;
    split(splInfo_sp, splInfo, is_any_of(" "));
    ul spl_num(splInfo_sp.size());
    vector<char> nuc_vec(spl_num);
    unordered_set<char> nuc_set;
    for (ul i = 0; i < spl_num; ++i) {
      vector<string> splInfo_sp_sp;
      split(splInfo_sp_sp, splInfo_sp.at(i), is_any_of(":"));
      string nuc(splInfo_sp_sp.at(0));
      string kACGT("ACGT");
      if (kACGT.find(nuc.at(0)) == string::npos || nuc.size() != 1) {
        nuc_vec.at(i) = 'n';
      } else {
        nuc_vec.at(i) = nuc.at(0);
        nuc_set.insert(nuc.at(0));
      }
    }
    if (nuc_set.size() >= 2) {
      cout << ctg_nom << '\t' << pos + 1;
      for (auto it = nuc_vec.begin(); it != nuc_vec.end(); ++it) {
        cout << '\t' << *it;
      }
      cout << '\n';
    }
  }
}
