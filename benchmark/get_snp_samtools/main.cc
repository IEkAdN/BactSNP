#include <iostream>
#include "two_matrix.h"

using std::cerr;

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "usage: 0 [.vcf_lst] [ref.fa]\n";
    return 1;
  } else {
    TwoMatrix::main(argv[1], argv[2]);
  }
  return 0;
}
