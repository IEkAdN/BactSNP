#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga(string PhenixFaNom, string RefFaNom) : kPhenixFaNom_(PhenixFaNom), kRefFaNom_(RefFaNom) {}
  ~Fuga() {}
  void Main();

 private:
  const string kPhenixFaNom_;
  const string kRefFaNom_;
  void ReadPhenixFa();
  void ReadRefFa();
  void Fa2Snp();
  map<string, string> Fa_;
  string RefChrId_;
};

#endif  // FUGA_H_
