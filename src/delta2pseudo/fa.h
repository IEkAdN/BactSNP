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

#ifndef FA_H
#define FA_H

#include "hoge.h"

class Fa {
 public:
  Fa() : kACGT("ACGT") {}
  ~Fa() {}
  const string kACGT;
  std::vector<std::string> id() const { return id_; }
  std::string seq(std::string id) const { return id_2_seq_.at(id); }
  std::unordered_map<std::string, std::string> id_2_seq() const { return id_2_seq_; }
  char nuc(std::string id, int64_t pos) const { return id_2_seq_.at(id).at(pos); }
  void ReadFa(const std::string& fa_nom);

 private:
  std::vector<std::string> id_;
  std::unordered_map<std::string, std::string> id_2_seq_;
};

#endif  // FA_H
