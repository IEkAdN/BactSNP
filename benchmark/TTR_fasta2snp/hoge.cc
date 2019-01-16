#include "hoge.h"

void split(vector<string>& s_sp, const string& s, const string& delim,
           bool token_compress) {
  s_sp.clear();  // same as boost::algorithm::split
  u64 pos_elemBeg(0);
  u64 pos_elemEnd(0);
  while ((pos_elemEnd = s.find_first_of(delim, pos_elemBeg))
         != std::string::npos) {
    if (token_compress) {
      // If delim characters are adjacent or
      // the first character of s is one of delim
      // (then, pos_elemEnd == pos_elemBeg),
      // empty string is not pushed back to s_sp.
      if (pos_elemEnd > pos_elemBeg) {
        s_sp.push_back(s.substr(pos_elemBeg, pos_elemEnd - pos_elemBeg));
      }
    } else {
      s_sp.push_back(s.substr(pos_elemBeg, pos_elemEnd - pos_elemBeg));
    }
    pos_elemBeg = pos_elemEnd + 1;
  }
  if (token_compress) {
    // If the last character of s is one of delim
    // (then, pos_elemBeg == s.length()),
    // empty string is not pushed back to s_sp.
    if (pos_elemBeg < s.length()) {
      s_sp.push_back(s.substr(pos_elemBeg));
    }
  } else {
    s_sp.push_back(s.substr(pos_elemBeg));
  }
}

void split_2(vector<string>& s_sp, const string& s, const string& delim,
             bool token_compress) {
  s_sp.clear();  // same as boost::algorithm::split
  u64 pos_elemBeg(0);
  u64 pos_elemEnd(0);
  while ((pos_elemEnd = s.find(delim, pos_elemBeg))
         != std::string::npos) {
    if (token_compress) {
      // If delim strings are adjacent or
      // the first part of s is delim
      // (then, pos_elemEnd == pos_elemBeg),
      // empty string is not pushed back to s_sp.
      if (pos_elemEnd > pos_elemBeg) {
        s_sp.push_back(s.substr(pos_elemBeg, pos_elemEnd - pos_elemBeg));
      }
    } else {
      s_sp.push_back(s.substr(pos_elemBeg, pos_elemEnd - pos_elemBeg));
    }
    pos_elemBeg = pos_elemEnd + delim.size();
  }
  if (token_compress) {
    // If the last part of s is delim (
    // (then, pos_elemBeg == s.length()),
    // empty string is not pushed back to s_sp.
    if (pos_elemBeg < s.length()) {
      s_sp.push_back(s.substr(pos_elemBeg));
    }
  } else {
    s_sp.push_back(s.substr(pos_elemBeg));
  }
}
