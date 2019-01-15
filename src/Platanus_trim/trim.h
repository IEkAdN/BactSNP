// -*- mode: c++ -*-
#ifndef __TRIM_H__
#define __TRIM_H__

#include "common.h"

class Trim : public BaseCommand {
public:

  static const int HIGH_QUALITY_REGION_TH;
  static const int Q_OPTION_LOWER_TH;

  typedef unordered_map<size_t, vector<int> > nlist;
  Trim();
  Trim(const unsigned l);
  virtual ~Trim() {}

  const char* usage() const;
  int exec();

  struct ADAPT_KMER : public unordered_map<unsigned,int> {
    static const int OFFSET=10000, len=11;

    ADAPT_KMER(string adp1="AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
               string adp2="CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT");
    ADAPT_KMER(string adp1);
  };

  struct INTERNAL_ADAPT : public map<vector<char>, int> {
    static const int len=18;
    int all_len;

    INTERNAL_ADAPT() {}
    INTERNAL_ADAPT(string adp, unsigned type);

  };

  struct Result
  {
    unsigned long long num_read_f;
    unsigned long long num_base_f;
    unsigned long long num_read_r;
    unsigned long long num_base_r;
    unsigned long long num_pair_or;
    unsigned long long num_pair_and;

    Result(): num_read_f(0), num_base_f(0), num_read_r(0), num_base_r(0), num_pair_or(0), num_pair_and(0) {}
    ~Result() = default;

    Result &operator+=(const Result &res)
    {
      num_read_f += res.num_read_f;
      num_base_f += res.num_base_f;
      num_read_r += res.num_read_r;
      num_base_r += res.num_base_r;
      num_pair_or += res.num_pair_or;
      num_pair_and += res.num_pair_and;
      return *this;
    }

    const Result operator+(const Result &res) const
    {
      Result ret;
      ret += res;
      return ret;
    }

  friend std::ostream &operator<<(std::ostream &os, const Result &res)
  {
    os << "NUM_OF_TRIMMED_READ(FORWARD) = " << res.num_read_f << "\n"
               << "NUM_OF_TRIMMED_BASE(FORWARD) = " << res.num_base_f << "\n"
               << "NUM_OF_TRIMMED_READ(REVERSE) = " << res.num_read_r << "\n"
               << "NUM_OF_TRIMMED_BASE(REVERSE) = " << res.num_base_r << "\n"
               << "NUM_OF_TRIMMED_PAIR(OR) = " << res.num_pair_or << "\n"
               << "NUM_OF_TRIMMED_PAIR(AND) = " << res.num_pair_and << "\n";
     return os;
  }


  };

    

  /**
   * @return number of trimmed adapter
   */
  
  template<class T, class Q=binder2nd<minus<char> > >
  inline void read_fastq_ecc(ifstream &ifs, vector<T> &v, nlist &n_list, Q qual_proc=bind2nd(minus<char>(),(char)64))
  {
    read_fastq_ecc(ifs, v, n_list, mem_fun_ref<void,vector<T>,const T&>(&vector<T>::push_back), qual_proc);
  };

  template<class A, class P, class Q>
  inline void read_fastq_ecc(ifstream &ifs, A &arg, nlist &n_list, P proc, Q qual_proc)
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
      check_n(line[0], n_list[num_of_reads-1]);
      transform(line[0].begin(), line[0].end(), seq.seq.begin(), char2bin4);
      transform(line[1].begin(), line[1].end(), seq.qual.begin(), qual_proc);
      proc(arg, seq);
    }
  };

  template<class T>
  inline void read_fasta_ecc(ifstream &ifs, vector<T> &v,nlist &n_list, char(*seq_fun)(char)=char2bin4)
  {
    read_fasta_ecc(ifs, v,n_list, mem_fun_ref<void,vector<T>,const T&>(&vector<T>::push_back), seq_fun);
  }

  template<class A, class P>
  inline void read_fasta_ecc(ifstream &ifs, A &arg,nlist &n_list, P proc, char(*seq_fun)(char)=char2bin4)
  {
    string  line[2];
    size_t  num_of_reads = 0;
    while (ifs && !ifs.eof() && num_of_reads++<MAX_PROCESSED_READS) {
      // read name and sequence
      SEQ  seq;
      if (!getline(ifs, seq.name) && line[0].empty()) {
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
      check_n(line[0], n_list[num_of_reads-1]);
      transform(line[0].begin(), line[0].end(), seq.seq.begin(), seq_fun);
      // call function and update name
      // seq.name = line[1];
      proc(arg, seq);
    if (!ifs.eof()) {  // rewind until name line
	ifs.putback('\n');
      foreach_r(line[1], c) {
        ifs.putback(*c);
      }
    }
    }
  };



protected:

  struct BLOCK {
    int qst, qed, sst, sed;
    char str, flag;
    BLOCK():qst(0),qed(0),sst(0),sed(0),str('/'),flag(-1) {}

    inline int operator- (BLOCK& rhs) { return (sst-rhs.sst-rhs.qst+qst); }
  };

  inline void check_n(string &line, vector<int> &n_list)
  {
    for(int i = 0; i < line.size(); i++) {
      if (line[i] == 'N') n_list.push_back(i);
    }
  }

  bool    _debug_mode;
  enum MODE {NOT, INTERNAL};
  MODE mode;
  int th_qv;

  int _exec(ifstream &kfs, ADAPT_KMER &adapt);
  int trim_adapter(ADAPT_KMER& adapt, SEQ &seq_f, SEQ &seq_r,vector<int> &n_list1,vector<int> &n_list2, Result &res, bool is_fastq=true);
  int _exec_internal(ifstream& kfs, ADAPT_KMER &adaptn, INTERNAL_ADAPT& adapt);
  int trim_internal(const INTERNAL_ADAPT& adapt, SEQ &seq_f, SEQ &seq_r, Result &res, bool is_fastq=true);
  int _trim_low_qual(SEQ &seq_f,SEQ &seq_r, vector<int> &n_list1,vector<int> &n_list2,Result &res);
  //template<class KMER> vector<char> _correct(KMER &kmer, vector<char> &seq, vector<unsigned> &error_pos);
  void _output(vector<SEQ> &seq1, vector<SEQ> &seq2, ofstream &ofsp, ofstream &ofss, nlist &n_list1, nlist &n_list2, int len,bool is_fastq);

};

#endif
