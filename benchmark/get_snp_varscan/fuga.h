#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga() {}
  ~Fuga() {}
  void Main(string result_nom, string isoLst_nom);

 private:
  void ReadResult(string f_nom);
  void ReadIsoLst(string f_nom);
  vector<string> isoLst_;
};

#endif  // FUGA_H_
