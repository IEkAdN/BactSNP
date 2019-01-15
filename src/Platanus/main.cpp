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
#include "assemble.h"
#include "scaffold.h"
#include "gapClose.h"

using std::cerr;
using std::endl;

//////////////////////////////////////////////////////////////////////////////////////
// show usage
//////////////////////////////////////////////////////////////////////////////////////
void Usage(void)
{
    cerr << "Usage: platanus Command [options]\n\nCommand: assemble, scaffold, gap_close" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////
// show version
//////////////////////////////////////////////////////////////////////////////////////
void ShowVersion(void)
{
    cerr << "Platanus version: " << platanus::ConstParam::VERSION << endl;
}

//////////////////////////////////////////////////////////////////////////////////////
// show input arguments
//////////////////////////////////////////////////////////////////////////////////////
void ShowArgs(const int argc, char **argv)
{
    for (int i = 0; i < argc; ++i) {
        cerr << argv[i] << " ";
    }
    cerr << endl << endl;
}



//////////////////////////////////////////////////////////////////////////////////////
// main function
//////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{
    // put_command_log(argc, argv, stderr);
	ShowVersion();
    ShowArgs(argc, argv);
    std::unique_ptr<BaseCommand> command;
    if (argc > 1) {
        if (strcmp(argv[1], "-v") == 0) {
            return 0;
        }
        if (strcmp(argv[1], "assemble") == 0) {
            command.reset(new Assemble());
        } else if (strcmp(argv[1], "scaffold") == 0) {
            command.reset(new Scaffold());
        } else if (strcmp(argv[1], "gap_close") == 0) {
            command.reset(new GapClose());
        }

        if (!command) {
            Usage();
            return 1;
        }
    } else {
        Usage();
        return 1;
    }
    if (!command->parseArgs(argc, argv)) {
        command->usage();
        return 1;
    }
    try {
        command->exec();
    } catch (std::string s) {
        cerr << s << endl;
    } catch (platanus::ErrorBase &e) {
        e.showErrorMessage();
        return e.getID();
    } catch (...) {
    }

  unsigned long long vmpeak = 0,vmhwm = 0;
  char proc_path[64];
  snprintf(proc_path, sizeof(proc_path), "/proc/%d/status", (int)getpid());    
  std::ifstream ifs(proc_path);
  std::string line;
  while (getline(ifs, line).good()) {
    const char* c = line.c_str();
    if (strlen(c) >= 6) {
      if (memcmp("VmPeak", c, 6) == 0) {
        vmpeak = atoi(c + 8);
      } else if (memcmp("VmHWM", c, 5)==0) {
        vmhwm = atoi(c + 8);
      }
    }
  }
  fprintf(stderr, "\n#### PROCESS INFORMATION ####\n");
  fprintf(stderr, "VmPeak:      %10.3f GByte\n", vmpeak/1048576.0);
  fprintf(stderr, "VmHWM:       %10.3f GByte\n", vmhwm/1048576.0);

  return 0;
}
