#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "usage: a [phynix.fasta] [ref.fasta]\n";
    return 1;
  } else {
    Fuga fuga(argv[1], argv[2]);
    fuga.Main();
  }
  return 0;
}
