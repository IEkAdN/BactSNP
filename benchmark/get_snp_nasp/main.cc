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
using std::stoul;
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
    vector<string> matrix_v;
    getline(matrix_f, matrix_l);
    split(matrix_v, matrix_l, is_any_of("\t"));
    unsigned iso_num(matrix_v.size() - 22);
    cout << "contig\tpos";
    for (unsigned i = 2; i <= iso_num + 1; ++i) {
      vector<string> matrix_v_v;
      split(matrix_v_v, matrix_v.at(i), is_any_of(":"));
      cout << "\t" << matrix_v_v.at(0);
    }
    cout << "\n";
    while (getline(matrix_f, matrix_l)) {
      split(matrix_v, matrix_l, is_any_of("\t"));
      set<char> alt_nuc;
      for (unsigned i = 2; i <= iso_num + 2; ++i) {
        char _alt_nuc(matrix_v.at(i).at(0));
        if (_alt_nuc == 'A' || _alt_nuc == 'C' ||
            _alt_nuc == 'G' || _alt_nuc == 'T') {
          alt_nuc.insert(_alt_nuc);
        }
      }
      if (alt_nuc.size() > 1) {
        vector<string> matrix_v_v;
        split(matrix_v_v, matrix_v.at(0), is_any_of(":"), boost::token_compress_on);
        cout << matrix_v_v.at(0) << "\t" << matrix_v_v.at(1);
        for (unsigned i = 2; i <= iso_num + 1; ++i) {
          cout << "\t" << matrix_v.at(i).at(0);
        }
        cout << "\n";
      }
    }
  }
  return 0;
}
