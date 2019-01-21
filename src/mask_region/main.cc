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

int main(int argc, const char* argv[]) {
  if (argc != 3) {
    cerr << "usage: a [in.fa] [in.masked.region]\n";
    return 1;
  } else {
    Fuga fuga;
    fuga.Main(argv[1], argv[2]);
  }
  return 0;
}
