#include <iostream>
#include "two_matrix.h"

using std::cerr;

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    cerr << "usage: 0 [.vcf]\n";
    return 1;
  } else {
    TwoMatrix two_matrix;
    two_matrix.main(argv[1]);
  }
  return 0;
}
