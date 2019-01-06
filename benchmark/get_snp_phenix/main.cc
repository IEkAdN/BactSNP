#include "fuga.h"

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    cerr << "usage: a []\n";
    return 1;
  } else {
    Fuga fuga(argv[1]);
    fuga.Main();
  }
  return 0;
}
