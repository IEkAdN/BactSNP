#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 7) {
    cerr << "usage: a [<prefix>_Internal_Nodes_True_alignment.FASTA of EvolveAGene] [<prefix>_True_alignment.FASTA] [1-to-1.delta] [ref_genome.fa] [root_genome.fa] [output directory]\n";
    return 1;
  } else {
    Fuga fuga(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6]);
    fuga.Main();
  }
  return 0;
}
