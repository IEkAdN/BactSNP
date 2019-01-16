#include "fuga.h"

void Fuga::Main() {
  ReadFa("Ref", kRefFaNom_);
  ReadIsoFaLst();
  SetRefGap();
  PrintSnp();
}

void Fuga::ReadIsoFaLst() {
  ifstream F(kIsoFaLstNom_);
  string L;
  while (getline(F, L)) {
    vector<string> LSp;
    split(LSp, L, "\t");
    ReadFa(LSp.at(0), LSp.at(1));
    Iso_.emplace_back(LSp.at(0));
  }
}

void Fuga::ReadFa(string Iso, string FNom) {
  Fa_[Iso];
  ifstream F(FNom);
  string L;
  while (getline(F, L)) {
    if (! L.empty()) {
      if (L.at(0) == '>') {
        if (Iso == "Ref") {
          ScafId_ = L.substr(1, L.find("_") - 1);
        }
        Fa_.at(Iso)[ScafId_] = "";
      } else {
        Fa_.at(Iso).at(ScafId_) += L;
      }
    }
  }
}

void Fuga::SetRefGap() {
  const string& RefScafSeq(Fa_.at("Ref").at(ScafId_));
  u32 ScafLen(RefScafSeq.size());
  RefGap_[ScafId_].resize(ScafLen, false);
  for (u32 j = 0; j != ScafLen; ++j) {
    if (RefScafSeq.at(j) == '-') {
      RefGap_.at(ScafId_).at(j) = true;
    }
  }
}

void Fuga::PrintSnp() {
  cout << "scaffold\tpos";
  for (auto i = Iso_.begin(); i != Iso_.end(); ++i) {
    cout << "\t" << *i;
  }
  cout << "\n";
  u32 ScafLen(Fa_.at("Ref").at(ScafId_).size());
  u32 WoGapPos(0);
  for (u32 j = 0; j < ScafLen; ++j) {
    if (! RefGap_.at(ScafId_).at(j)) {
      unordered_set<char> Nuc;
      for (auto k = Iso_.begin(); k != Iso_.end(); ++k) {
        char _Nuc(Fa_.at(*k).at(ScafId_).at(j));
        if (kACGT_.find(_Nuc) != string::npos) {
          Nuc.insert(_Nuc);
        }
      }
      if (Nuc.size() >= 2) {
        cout << ScafId_ << "\t" << WoGapPos + 1;
        for (auto k = Iso_.begin(); k != Iso_.end(); ++k) {
          cout << "\t" << Fa_.at(*k).at(ScafId_).at(j);
        }
        cout << "\n";
      }
      ++WoGapPos;
    }
  }
}
