#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "usage: a [base_genome_aligned.fa] [simulated_genome_aligned.fasta_list]\n";
    return 1;
  } else {
    Fuga fuga(argv[1], argv[2]);
    fuga.Main();
  }
  return 0;
}
