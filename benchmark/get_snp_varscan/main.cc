#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "usage: a [in.VarScan.result] [in.isolate_list]\n";
    return 1;
  } else {
    Fuga fuga;
    fuga.Main(argv[1], argv[2]);
  }
  return 0;
}
