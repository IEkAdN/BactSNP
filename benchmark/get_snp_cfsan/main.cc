#include <iostream>
#include "hoge.h"

using std::cerr;
using std::stod;
using std::stoul;
using std::stoull;

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    cerr << "usage: a [in.vcf]\n";
    return 1;
  } else {
    Hoge hoge;
    hoge.Main(argv[1]);
  }
  return 0;
}
