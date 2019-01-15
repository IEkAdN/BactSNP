// -*- mode: c++ -*-
#ifndef __COMMON_H__
#define __COMMON_H__

#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <functional>
#include <string.h>
#include <omp.h>
#include <getopt.h>
#include "logger.h"
using namespace std;



#define foreach(CNTR,I) for(auto I=(CNTR).begin(),I_e=(CNTR).end(); I!=I_e; I++)
#define foreach_r(CNTR,I) for(auto I=(CNTR).rbegin(),I_e=(CNTR).rend(); I!=I_e; I++)

// char2bin5: A=>0, C=>1, G=>2, T=>3, N=>4
// char2bin4: A=>0, C=>1, G=>2, T=>3, N=>0
// ulong2bin: 0x1=>0, 0x2=>1, 0x4=>2, 0x8=>3
// bin2char:  0=>A, 1=>C, 2=>G, 3=>T, 4=>N
inline char char2bin5(char c) { return ".\x0.\x1\x3..\x2......\x4"[c&0xF]; }
inline char char2bin4(char c) { return ".\x0.\x1\x3..\x2......\x0"[c&0xF]; }
/*
inline char char2bin4(char c) {
  return
    "................................................................."
    "\x0.\x1...\x2......\x0.....\x3............"
    "\x0.\x1...\x2......\x0.....\x3"[c];
}
*/
#define ulong2bin(u) ".\x0\x1.\x2...\x3"[(u)&0xF]
#define bin2char(b)  "ACGTN"[b]
struct FASTQ_PARAM {
    static const int FASTQ_OLD;
    static const int FASTQ_NEW;
};


/** insert comma into integer */
string number_format(size_t num);


class BaseCommand {
public:
  enum FILETYPE { UNKNOWN, FASTA, FASTQ33, FASTQ64 };


  const char*                     short_options;
  unordered_map<char,string>      exec_option;

  static vector<pair<FILETYPE,string> >  input_file;
  static Logger                          logger;
  static int                             NUM_OF_THREADS;
  static size_t                          MAX_PROCESSED_READS;

  struct SEQ {
    vector<char>  seq, qual;
    string        name;

    void put(ostream& os, vector<int> &n_list, bool with_quality=false) {
      os << name.c_str() << '\n';
      int i = 0;
      vector<int>::iterator it;
      if (!n_list.empty())
	    for(it=n_list.begin();(it!=n_list.end())&&(*it<0);it++);
      foreach (seq, s) {
        if(!n_list.empty() && it != n_list.end() && i == *it) {
          os.put('N');
          ++it;
        } else
          os.put(bin2char(*s));
        i++;
      }
      os << '\n';
      if (with_quality) {
        os << '+' << name.substr(1).c_str() << '\n';
        foreach (qual, s) {
          os.put((char)(*s)+33);
        }
        os << '\n';
      }
    }
  };

  BaseCommand() {}
  virtual ~BaseCommand() {}

  virtual bool parse_args(int argc, char* argv[]);
  virtual const char* usage() const=0;
  virtual int exec()=0;


  template<class T, class Q=binder2nd<minus<char> > >
  inline void read_fastq(ifstream &ifs, vector<T> &v, Q qual_proc=bind2nd(minus<char>(),(char)64))
  {
    read_fastq(ifs, v, mem_fun_ref<void,vector<T>,const T&>(&vector<T>::push_back), qual_proc);
  }

  template<class A, class P, class Q>
  inline void read_fastq(ifstream &ifs, A &arg, P proc, Q qual_proc)
  {
    string  line[2];
    size_t  num_of_reads = 0;
    while (ifs && !ifs.eof() && num_of_reads++<MAX_PROCESSED_READS) {
      // read sequece and quality
      SEQ  seq;
      if (!getline(ifs, seq.name)       // name
          || !getline(ifs, line[0])     // sequence (=read)
          || !getline(ifs, line[1])     // name (skipped)
          || !getline(ifs, line[1])) {  // quality
        break;
      }
      // transform binary from charactor 'ACGTN' and call member function
      seq.seq.resize(line[0].size());
      seq.qual.resize(line[1].size());
      transform(line[0].begin(), line[0].end(), seq.seq.begin(), char2bin4);
      transform(line[1].begin(), line[1].end(), seq.qual.begin(), qual_proc);
      proc(arg, seq);
    }
  };

  template<class T>
  inline void read_fasta(ifstream &ifs, vector<T> &v, char(*seq_fun)(char)=char2bin4)
  {
    read_fasta(ifs, v, mem_fun_ref<void,vector<T>,const T&>(&vector<T>::push_back), seq_fun);
  }

  template<class A, class P>
  inline void read_fasta(ifstream &ifs, A &arg, P proc, char(*seq_fun)(char)=char2bin4)
  {
    string  line[2];
    size_t  num_of_reads = 0;
    while (ifs && !ifs.eof() && num_of_reads++<MAX_PROCESSED_READS) {
      // read name and sequence
      SEQ  seq;
      if (line[0].empty() && !getline(ifs, seq.name)) {
        break;
      }
      line[0].clear();
      while (getline(ifs, line[1])) {  // read sequence
        if (line[1].empty()) {  // empty line is skipped
          continue;
        }
        if (!strchr("ACGTN ", line[1][0]) || ifs.eof()) {
          break;   // this line is name
        }
        line[0] += line[1];
      }
      if (line[0].empty()) {
        break;
      }
      // remove space and transform binary from charactor 'ACGTN'
      std::remove(line[0].begin(), line[0].end(), ' ');
      seq.seq.resize(line[0].size());
      transform(line[0].begin(), line[0].end(), seq.seq.begin(), seq_fun);
      // call function and update name
      seq.name = line[1];
      proc(arg, seq);
    }
    if (!ifs.eof()) {  // rewind until name line
      ifs.putback('\n');
      foreach_r(line[1], c) {
        ifs.putback(*c);
      }
    }
  };

protected:
  bool _check_fastq_format(ifstream& ifs);
  FILETYPE _guess_filetype(string& path);
};

#endif
