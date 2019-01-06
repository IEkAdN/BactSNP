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
  Hoge() : pos_(0), gt_idx_(0), cov_idx_(0), iso_num_(0) {}
  ~Hoge() {}
  void Main(const std::string& f_nom);

 private:
  void ReadVcf(const std::string& f_nom);
  void Filter(const std::vector<std::string>& l_sp);
  void ReadFormatVal(const std::vector<std::string>& l_sp,
                     unsigned long iso_idx);
  unsigned long GetGT(const std::vector<std::string>& formatVal_sp);
  void SetIsoNuc(unsigned long gt);
  void ResetVar();
  std::string cntg_id_;
  unsigned long pos_;
  char ref_nuc_;
  char alt_nuc_;
  std::string filter_;
  std::string info_;
  std::vector<std::string> info_sp_;
  std::string svtype_;
  std::string format_;
  std::vector<std::string> format_sp_;
  unsigned long gt_idx_;
  unsigned long cov_idx_;
  unsigned long iso_num_;
  std::string formatVal_;
  std::unordered_set<unsigned long> gt_set_;
  std::vector<char> iso_nuc_;
  bool is_there_cov0AlleleSpl_;
};

#endif  // HOGE_H_
