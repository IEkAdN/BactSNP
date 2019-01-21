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

#include <iostream>

#include "base_caller.h"

using std::cerr;
using std::stoi;
using std::stoll;

int main(int argc, const char* argv[]) {
  if (argc != 5) {
    cerr << "usage: 0 [in.ref.fa] [in.qry.fa] [in.delta] "
         << "<length near indel to be masked>\n";
    return 1;
  } else {
    BaseCaller base_caller(stoll(argv[4]));
    base_caller.Main(argv[1], argv[2], argv[3]);
  }
  return 0;
}
