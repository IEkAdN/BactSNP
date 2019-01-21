/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "fuga.h"

void Fuga::Main(string sam_nom, string R1_fq_nom, string R2_fq_nom) {
  ReadSam(sam_nom, R1_fq_nom, R2_fq_nom);
}

void Fuga::ReadSam(string sam_nom, string R1_fq_nom, string R2_fq_nom) {
  ifstream sam(sam_nom);
  ofstream R1_fq(R1_fq_nom);
  ofstream R2_fq(R2_fq_nom);
  string sam_l;
  // true -> R1, false -> R2
  bool read_1_or_2(false);
  while (getline(sam, sam_l)) {
    if (sam_l.at(0) != '@') {
      read_1_or_2 = ! read_1_or_2;
      vector<string> sam_l_sp;
      split(sam_l_sp, sam_l, "\t");
      string id(sam_l_sp.at(0));
      ul flag(stoul(sam_l_sp.at(1)));
      string seq(sam_l_sp.at(9));
      if (flag & 0x10) {
        RevComp(&seq);
      }
      if (read_1_or_2) {
        R1_fq << "@" << id << "/1\n"
               << seq << "\n"
               << "+\n"
               << string(seq.size(), 'I') << "\n";
      } else {
        R2_fq << "@" << id << "/2\n"
               << seq <<"\n"
               << "+\n"
               << string(seq.size(), 'I') << "\n";
      }
    }
  }
}

void Fuga::RevComp(string* seq) {
  string rev_comp_seq("");
  rev_comp_seq.reserve(seq->size());
  for (auto it = seq->rbegin(); it != seq->rend(); ++it) {
    if (*it == 'A' || *it == 'a') {
      rev_comp_seq += 'T';
    } else if (*it == 'C' || *it == 'c') {
      rev_comp_seq += 'G';
    } else if (*it == 'G' || *it == 'g') {
      rev_comp_seq += 'C';
    } else if (*it == 'T' || *it == 't') {
      rev_comp_seq += 'A';
    } else {
      rev_comp_seq += 'N';
    }
  }
  *seq = rev_comp_seq;
}
