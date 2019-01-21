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

#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga {
 public:
  Fuga() : kValidNuc("ACGTacgt") {}
  ~Fuga() {}
  void Main(const std::string& f_nom);

 private:
  void ReadFaLst(const std::string& f_nom);
  void SetCntgInfo(const std::string& f_nom);
  void ReadFa(const std::string& iso_nom, const std::string& f_nom);
  void Print();
  const string kValidNuc;
  std::vector<std::string> iso;
  std::unordered_map<std::string, unsigned long> cntg_info;
  std::unordered_map<std::string,
                     std::unordered_map<std::string, std::string>*> cntg_seq;
};

#endif  // FUGA_H_
