#ifndef HOGE_H_
#define HOGE_H_

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class Hoge {
 public:
  Hoge() {}
  ~Hoge() {}
  void Main(const std::string& f_nom);

 private:
  void ReadVcf(const std::string& f_nom);
  bool LineIsExpectedFmt(const std::string& l,
                         const std::vector<std::string>& l_sp);
  void PrintInterIsoSNP(const std::string& l,
                        const std::vector<std::string>& l_sp);
  unsigned long GetGT(const std::vector<std::string>& l_sp,
                      unsigned long gt_idx, unsigned long iso_idx);
  void SetIsoNuc(unsigned long gt);
};

#endif  // HOGE_H_
