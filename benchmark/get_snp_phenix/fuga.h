#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga(string FNom) : kFNom_(FNom) {}
  ~Fuga() {}
  void Main();

 private:
  const string kFNom_;
  void Read();
  void Fa2Snp();
  map<string, string> Fa_;
};

#endif  // FUGA_H_
