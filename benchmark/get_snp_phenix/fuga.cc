#include "fuga.h"

void Fuga::Main() {
  Read();
  Fa2Snp();
}

void Fuga::Read() {
  ifstream F(kFNom_);
  string L;
  getline(F, L);
  getline(F, L);
  string FaId;
  while (getline(F, L)) {
    if (! L.empty()) {
      if (L.at(0) == '>') {
        FaId = L.substr(1, L.find(" ") - 1);
        Fa_[FaId] = "";
      } else {
        Fa_.at(FaId) += L;
      }
    }
  }
}

void Fuga::Fa2Snp() {
  cout << "pos";
  for (auto i = Fa_.begin(); i != Fa_.end(); ++i) {
    cout << "\t" << i->first;
  }
  cout << "\n";
  u32 SeqLen(Fa_.begin()->second.size());
  for (u32 i = 0; i < SeqLen; ++i) {
    unordered_set<char> Nuc;
    for (auto j = Fa_.begin(); j != Fa_.end(); ++j) {
      Nuc.insert(j->second.at(i));
    }
    // output of vcf2fasta contains only ACGTN-
    Nuc.erase('-');
    Nuc.erase('N');
    if (Nuc.size() >= 2) {
      cout << i + 1;
      for (auto j = Fa_.begin(); j != Fa_.end(); ++j) {
        cout << "\t" << j->second.at(i);
      }
      cout << "\n";
    }
  }
}
