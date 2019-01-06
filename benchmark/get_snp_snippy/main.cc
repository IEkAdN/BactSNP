#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using std::cerr;
using std::cout;
using std::ifstream;
using std::set;
using std::string;
using std::vector;
using boost::algorithm::is_any_of;
using boost::algorithm::split;

int main(int argc, const char* argv[]) {
  if (argc != 2) {
    cerr << "usage: 0 [.matrix]\n";
    return 1;
  } else {
    ifstream matrix_f(argv[1]);
    string matrix_l = "";
    getline(matrix_f, matrix_l);
    cout << matrix_l << '\n';
    while (getline(matrix_f, matrix_l)) {
      vector<string> matrix_v;
      split(matrix_v, matrix_l, is_any_of("\t"));
      set<char> alt_nuc;
      for (unsigned i = 2; i < matrix_v.size(); ++i) {
        char _alt_nuc(matrix_v.at(i).at(0));
        if (_alt_nuc == 'A' || _alt_nuc == 'C' ||
            _alt_nuc == 'G' || _alt_nuc == 'T') {
          alt_nuc.insert(_alt_nuc);
        }
      }
      if (alt_nuc.size() > 1) cout << matrix_l << '\n';
    }
  }
  return 0;
}
