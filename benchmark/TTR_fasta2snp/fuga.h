#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga(string RefFaNom, string IsoFaLstNom)
    : kRefFaNom_(RefFaNom), kIsoFaLstNom_(IsoFaLstNom), kACGT_("ACGT") {}
  ~Fuga() {}
  void Main();
  const string kRefFaNom_;
  const string kIsoFaLstNom_;
  const string kACGT_;

 private:
  void ReadFa(string Iso, string FNom);
  void ReadIsoFaLst();
  void SetRefGap();
  void PrintSnp();
  unordered_map<string, string> RefFa_;
  string ScafId_;
  // isolate -> scaf ID -> seq
  unordered_map<string, unordered_map<string, string> > Fa_;
  unordered_map<string, vector<bool> > RefGap_;
  vector<string> Iso_;
};

#endif  // FUGA_H_
