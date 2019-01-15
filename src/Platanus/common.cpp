/*
Copyright (C) 2014 Itoh Laboratory, Tokyo Institute of Technology

This file is part of Platanus.

Platanus is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

Platanus is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with Platanus; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "common.h"

//////////////////////////////////////////////////////////////////////////////////////
// static constant
//////////////////////////////////////////////////////////////////////////////////////
const unsigned long long platanus::ConstParam::MAX_READ_LEN = 5000;
const unsigned platanus::ConstParam::SCAFFOLD_HASH_OVERLAP = 32;
const unsigned platanus::ConstParam::OUTPUT_LINE_LENGTH = 80;
const unsigned platanus::ConstParam::MAX_FILE_NUM = 100;
const unsigned platanus::ConstParam::MAX_FILE_LEN = 200;
const unsigned platanus::ConstParam::MAX_THREAD = 100;
const std::string platanus::ConstParam::VERSION = "1.2.4";
const double platanus::ConstParam::DOUBLE_HASH_MAX_LOAD_FACTOR = 0.9;
std::string platanus::globalTmpFileDir = ".";



