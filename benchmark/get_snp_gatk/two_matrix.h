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
  void ReadVcf(const std::string& vcf_na);
  void main(const std::string& vcf_na);

 private:
  std::vector<std::string> iso_id_;
  DISALLOW_COPY_AND_ASSIGN(TwoMatrix);
};

#endif  // VCF_2MATRIX_02_TWO_MATRIX_H_
