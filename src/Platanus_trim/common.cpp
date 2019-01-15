
#include "common.h"




// initialize static members
vector<pair<BaseCommand::FILETYPE,string> >  BaseCommand::input_file;
Logger  BaseCommand::logger;
int     BaseCommand::NUM_OF_THREADS      = 0;
size_t  BaseCommand::MAX_PROCESSED_READS = 1;  // 1 = 100000
const int FASTQ_PARAM::FASTQ_OLD = 64;
const int FASTQ_PARAM::FASTQ_NEW = 33;



/** insert comma into integer */
string number_format(size_t num)
{
  if (num==0) {
    return string("0");
  }
  string buf;
  for (int i=0; num>0; num/=10,i++) {
    buf += '0'+(num%10);
    if (i%3==2 && num>9) {
      buf += ',';
    }
  }
  char c;
  string::iterator b=buf.begin(),e=buf.end();
  while ((b!=e)&&(b!=--e)) {
    c = *b;
    *b++ = *e;
    *e = c;
  }
  return buf;
}


bool BaseCommand::parse_args(int argc, char* argv[])
{
  // read arguments.
  optind = 1;
  while (1) {
    int c = getopt(argc, argv, short_options);
    if (c<0) {
      break;
    }
    if (c==0 || c=='?') {
      // logger.error() << usage() << endl;
      return false;
    }
    exec_option[(char)c] = (optarg ? optarg : "true");
  }

  // set input file pathes.
  string    path;
  FILETYPE  ftype;
  if (exec_option['i'] != "") {
    ifstream ifs(exec_option['i'].c_str());
    if (!ifs) {
      logger.error() << "Wrong file name: " << exec_option['i'] << endl;
      return false;
    }
    while (getline(ifs, path)) {
      ftype = _guess_filetype(path);
      if (ftype == FILETYPE::UNKNOWN) {
        logger.error() << "Invalid format or cannot open: " << path << endl;
        return false;
      }
      input_file.push_back(make_pair(ftype,path));
    }
    if (input_file.empty()) {
      logger.notice() << "None is valid file path: " << path << endl;
    }
  }
  while (optind < argc) {
    path = argv[optind++];
    ftype = _guess_filetype(path);
    if (ftype == FILETYPE::UNKNOWN) {
      logger.error() << "Invalid format or cannot open: " << path << endl;
      return false;
    }
    input_file.push_back(make_pair(ftype,path));
  }
  if (input_file.empty()) {
    return false;
  }

  // set number of threads
  if (exec_option.find('t')==exec_option.end()) {
    NUM_OF_THREADS = 1;
  } else {
    NUM_OF_THREADS = atoi(exec_option['t'].c_str());
    if (NUM_OF_THREADS<=0) {
      logger.error() << "Number of threads must be positive: " << NUM_OF_THREADS << endl;
      return false;
    }
  }

  // set max processed reads
  if (exec_option.find('n')==exec_option.end()) {
    MAX_PROCESSED_READS = 1;
  } else {
    MAX_PROCESSED_READS = atoi(exec_option['n'].c_str());
    if (MAX_PROCESSED_READS<=0) {
      logger.error() << "Max processed reads must be positive: " << MAX_PROCESSED_READS << endl;
      return false;
    }
  }
  MAX_PROCESSED_READS *= 100000;

  return true;
}

bool BaseCommand::_check_fastq_format(ifstream& ifs)
{
  int i;
  string line[4];
  while (ifs && !ifs.eof() && getline(ifs, line[0])
    && getline(ifs, line[1])
    && getline(ifs, line[2])
    && getline(ifs, line[3])) {
    for (i = 0; i < line[3].size(); ++i) {
      if ((int)line[3][i] - 74 > 0) { // check FASTQ ASCII -64
        return false;
      } else if ((int)line[3][i] - 59 < 0) { // check FASTQ ASCII -33
        return true;
      }
    }
  }
  return true;
}




BaseCommand::FILETYPE BaseCommand::_guess_filetype(string& path)
{
  string line[4];
  ifstream ifs(path.c_str());
  if (ifs
      && getline(ifs, line[0])
      && getline(ifs, line[1])
      && getline(ifs, line[2])
      && getline(ifs, line[3])) {

    if (line[0][0]=='>'
        && strspn(line[1].c_str(), "ACGTN")==line[1].size()) {
      return FILETYPE::FASTA;
    }
    if (line[0][0]=='@'
        && strspn(line[1].c_str(), "ACGTN")==line[1].size()
        && line[2][0]=='+') {
        return _check_fastq_format(ifs) ? FILETYPE::FASTQ33 : FILETYPE::FASTQ64;
    }
  }
  return FILETYPE::UNKNOWN;
}
