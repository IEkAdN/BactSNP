#ifndef VCF_2MATRIX_02_TWO_MATRIX_H_
#define VCF_2MATRIX_02_TWO_MATRIX_H_

#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
  TypeName(const TypeName&);               \
  void operator=(const TypeName&)

class TwoMatrix {
 public:
  TwoMatrix() {}
  ~TwoMatrix() {}
  static void ReadRef(const std::string& ref_na);
  static void ReadVcfLst(const std::string& vcf_lst_na);
  static void ReadVcf(const std::string& iso_id, const std::string& vcf_na);
  static void PrintMatrix();
  static void PrintMatrixLine(const std::string& contig, int64_t pos,
                              char ref_nuc);
  static void main(const std::string& vcf_lst_na, const std::string& ref_na);

 private:
  static std::vector<std::string> iso_id_;
  static std::map<std::string, std::string> contig_2_ref_seq;
  static std::map<std::string, std::map<std::string, std::string>>
      iso_id_2_contig_2_seq;
  DISALLOW_COPY_AND_ASSIGN(TwoMatrix);
};

#endif  // VCF_2MATRIX_02_TWO_MATRIX_H_
