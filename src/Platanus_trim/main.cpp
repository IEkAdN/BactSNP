
#include <unistd.h>
#include <sys/times.h>
#include <sys/types.h>
#include "trim.h"


int main (int argc, char* argv[])
{
 ios::sync_with_stdio(false);
  //
  // select target command
  //
  BaseCommand *cmd = NULL;
  if (RUN_MODE == "internal") {
      cmd = new Trim(static_cast<unsigned>(1));
    }
    else {
      cmd = new Trim();
    }
  if (!cmd) {
    cerr << "Fatally faied!!" << endl << endl
         << "Please re-compile again!!" << endl << endl;
    return 1;
  }

  //
  // parse arguments
  //
  if (!cmd->parse_args(argc, argv)) {
    cerr << cmd->usage() << endl;
    delete cmd;
    return 1;
  }

  //
  // execute command
  //
  time_t     start_at, end_at;
  struct tms t;
  time(&start_at);
  try {
    cmd->exec();
  } catch(...) {
    fprintf(stderr, "catch exception...\n");
  }
  delete cmd;
  time(&end_at);
  times(&t);

  //
  // output process informations
  //
  unsigned long vmpeak=0,vmhwm=0;
  char proc_path[64];
  snprintf(proc_path, sizeof(proc_path), "/proc/%d/status", (int)getpid());    
  ifstream ifs(proc_path);
  string line;
  while (getline(ifs, line).good()) {
    const char* c = line.c_str();
    if (strlen(c)>=6) {
      if (memcmp("VmPeak",c,6)==0) {
        vmpeak = atoi(c+8);
      } else if (memcmp("VmHWM",c,5)==0) {
        vmhwm = atoi(c+8);
      }
    }
  }
  fprintf(stderr, "\n#### PROCESS INFORMATION ####\n");
  fprintf(stderr, "User Time:     %10.2f min\n", t.tms_utime/6000.0);
  fprintf(stderr, "System Time:   %10.2f min\n", t.tms_stime/6000.0);
  fprintf(stderr, "VmPeak:      %10.3f GByte\n", vmpeak/1048576.0);
  fprintf(stderr, "VmHWM:       %10.3f GByte\n", vmhwm/1048576.0);
  fprintf(stderr, "Execution time:%10.2f min\n\n", difftime(end_at, start_at)/60.0);

  return 0;
}
