
#include "trim.h"
using namespace std;

const int Trim::HIGH_QUALITY_REGION_TH = 11;
const int Trim::Q_OPTION_LOWER_TH = 0;


Trim::Trim()
{
    short_options = "i:q:t:1:2:l:f";
    exec_option['i'] = "";
    exec_option['q'] = "15";
    exec_option['t'] = "1";
    exec_option['1'] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    exec_option['2'] = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT";
    exec_option['l'] = "25";
    exec_option['f'] = "false";
    mode = NOT;
}

Trim::Trim(const unsigned l)
{
    short_options = "i:q:t:1:2:a:b:rl:f";
    exec_option['q'] = "15";
    exec_option['i'] = "";
    exec_option['t'] = "1";
    exec_option['1'] = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
    exec_option['2'] = "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT";
    exec_option['a'] = "2";
    exec_option['b'] = "";
    exec_option['r'] = "false";
    exec_option['l'] = "25";
    exec_option['f'] = "false";
    mode = INTERNAL;
}



const char* Trim::usage() const
{
    switch (mode) {
	  case NOT: {
			    return "platanus_trim version 1.0.7\n\n" \
				  "usage: platanus_trim [options]\n\n"          \
				  "options:\n"                                                        \
				  "  -i str    List of input files (default NULL)\n"                  \
				  "  ------\n"                                                        \
				  "  -q int    Quality cutoff value (default 15)\n"                   \
				  "  -l int    Output length cutoff value  (default 25)\n"                   \
				  "  -f        Not remove read even if the pair is too short   \n"                   \
				  "  -t int    Max number of threads (default 1)\n"                   \
				  "  -1 str    adaptor 1 (default AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)\n"                                \
				  "  -2 str    adaptor 2 (default CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT)\n"                                \
				  "";
			    break;
			}
	  case INTERNAL: {
				   return "platanus_internal_trim version 1.0.7\n\n" \
					 "usage: platanus_trim [options]\n\n"          \
					 "options:\n"                                                        \
					 "  -i str    List of input files (default NULL)\n"                  \
					 "  ------\n"                                                        \
					 "  -t int    Max number of threads (default 1)\n"                   \
					 "  -a int    Select internal adaptor type (default 2)\n"                                 \
					 "            1: 454 adaptor (ATCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACGT)\n" \
					 "            2: Nextera adaptor (CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG)\n"    \
					 "            3: Solid adaptor (TGCTGTACCGTACATCCGCCTTGGCCGTACAGCAG)\n"         \
					 "  -b str    Set internal adaptor seq (if you can use original internal adaptor without -a option)\n"                                \
					 "  -q int    Quality cutoff value (default 15)\n"                   \
					 "  -l int    Output length cutoff value  (default 25)\n"                   \
					 "  -f        Not remove read even if the pair is too short   \n"                   \
					 "  -1 str    adaptor 1 (default AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)\n"                                \
					 "  -2 str    adaptor 2 (default CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT)\n"                                \
					 "";
				   break;
			     }
    }
}


Trim::ADAPT_KMER::ADAPT_KMER(string adp1, string adp2)
{
    unsigned fwd, rev, mask=(0x1ull<<(len*2))-0x1ull;
    int i, j, size, offset=0;
    for (string adp=adp1; offset<=this->OFFSET; offset+=this->OFFSET,adp=adp2) {
	  size = adp.size();
	  for (i=0; i<size-len+1; i++) {
		for (j=0; j<len; j++) {
		    fwd = (fwd<<2) | (char2bin4(adp[i+j]));
		    rev = (rev<<2) | (char2bin4(adp[i+len-1-j])^0x3);
		}
		(*this)[fwd & mask] = size - i - len + 1 + offset;
		(*this)[rev & mask] = i - size + len - 1 + offset;
	  }
    }
}

Trim::INTERNAL_ADAPT::INTERNAL_ADAPT(string adp, unsigned type)
{
    string adapt1 = "ATCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACGT";
    string adapt2 = "CTGTCTCTTATACACATCTAGATGTGTATAAGAGACAG";
    string adapt3 = "TGCTGTACCGTACATCCGCCTTGGCCGTACAGCAG";

    if (type == 1) {
	  adp = adapt1;
    } else if (type == 2) {
	  adp = adapt2;
    } else if (type == 3) {
	  adp = adapt3;
    } else if (adp.length() == 0) {
	  adp = adapt2;
    }
    all_len = adp.length();
    vector<char> fwd, rev;
    for (int i = 0; i < len; ++i) {
	  fwd.push_back(char2bin4(adp[i]));
	  rev.push_back(char2bin4(adp[adp.length() - 1 - i])^0x3);
    }
    for (int j = 0; j < len; ++j) {
	  vector<char> fwd_t = fwd;
	  vector<char> rev_t = rev;
	  for (char base = 0; base < 4; ++base) {
		fwd_t[j] = base;
		rev_t[j] = base;
		/*            for (auto it = fwd_t.begin(), end =fwd_t.end(); it != end; ++it) {
				  cout << (bin2char(*it));
				  } cout << endl;
				  for (auto it = rev_t.begin(), end =rev_t.end(); it != end; ++it) {
				  cout << (bin2char(*it));
				  } cout << endl;*/
		(*this)[fwd_t] = 1;
		(*this)[rev_t] = -1;
	  }
    }

    for (int i = len; i < adp.length(); ++i) {
	  fwd.push_back(char2bin4(adp[i]));
	  rev.push_back(char2bin4(adp[adp.length() - 1 - i])^0x3);
	  for (int j = 0; j < fwd.size(); ++j) {
		vector<char> fwd_t = fwd;
		vector<char> rev_t = rev;
		for (char base = 0; base < 4; ++base) {
		    fwd_t[j] = base;
		    rev_t[j] = base;
		    /*            for (auto it = fwd_t.begin(), end =fwd_t.end(); it != end; ++it) {
					cout << (bin2char(*it));
					} cout << endl;
					for (auto it = rev_t.begin(), end =rev_t.end(); it != end; ++it) {
					cout << (bin2char(*it));
					} cout << endl;
					*/
		    (*this)[fwd_t] = 1;
		    (*this)[rev_t] = -1;
		}
	  }
    }
    for (int i = len; i < adp.length(); ++i) {
	  fwd.erase(fwd.begin());
	  rev.erase(rev.begin());
	  for (int j = 0; j < fwd.size(); ++j) {
		vector<char> fwd_t = fwd;
		vector<char> rev_t = rev;
		for (char base = 0; base < 4; ++base) {
		    fwd_t[j] = base;
		    rev_t[j] = base;
		    (*this)[fwd_t] = 1;
		    (*this)[rev_t] = -1;
		}
	  }
    }

}




int Trim::exec()
{
    // load kmer information
    int          k, bucket;
    double       cov;
    ifstream     ifs;
    ADAPT_KMER   adapt;
    INTERNAL_ADAPT internal_adapt;
    th_qv = atoi(exec_option['q'].c_str());
    if (th_qv < Q_OPTION_LOWER_TH) {
	  logger.error() << "ERROR: -q option value must be >= " << Q_OPTION_LOWER_TH << "!!" << endl;
	  return 1;
    }
    switch (mode) {
	  case NOT: {
			    adapt = ADAPT_KMER(exec_option['1'], exec_option['2']);
			    // init parameters
			    _debug_mode = false;
			    logger.notice() << "Running with trim adapter mode";
			    logger.notice() << endl;

			    // execute
			    if (input_file.size()%2 == 1) {
				  logger.error() << "number of input files must be even." << endl;
			    } else {
				  return _exec(ifs, adapt);
			    }
			    break;
			}
	  case INTERNAL: {
				   adapt = ADAPT_KMER(exec_option['1'], exec_option['2']);
				   internal_adapt = INTERNAL_ADAPT(exec_option['b'], atoi(exec_option['a'].c_str()));
				   // init parameters
				   _debug_mode = false;
				   logger.notice() << "Running with trim internal adapter mode";
				   logger.notice() << endl;

				   // execute
				   if (input_file.size()%2 == 1) {
					 logger.error() << "number of input files must be even." << endl;
				   } else {
					 return _exec_internal(ifs, adapt, internal_adapt);
				   }

				   break;
			     }
	  default: {
			   throw;
			   break;
		     }
    }

    return 0;
}


int Trim::_exec(ifstream& kfs, Trim::ADAPT_KMER& adapt)
{
    omp_set_num_threads(NUM_OF_THREADS);

    vector<Result> adapt_result(NUM_OF_THREADS);
    vector<Result> qual_result(NUM_OF_THREADS);
    bool is_fastq;
    int output_th_len=atoi(exec_option['l'].c_str());
    // load k-mers and put partial on memory.


    for (auto f=input_file.begin(),f_e=input_file.end(); f!=f_e; f+=2) {
	  // generate output file pathes
	  ofstream dbg;
	  string fn[4]={(*f).second, (*(f+1)).second};
	  fn[2] = fn[0] + ".trimmed";
	  fn[3] = fn[1] + ".trimmed";

	  // open files
	  int i;
	  ifstream ifs[2];
	  ofstream ofs[2];
	  const int bufsize=1000000;
	  char buf[bufsize];
	  for (i=0; i<2; i++) {
		ifs[i].open(fn[i].c_str());
		ofs[i].open(fn[i+2].c_str());
		if (!ifs[i].good()) {
		    logger.error() << "ERROR: cannot open: " << fn[i] << endl;
		}
		if (!ofs[i].good()) {
		    logger.error() << "ERROR: cannot open: " << fn[i+2] << endl;
		}
		ifs[i].rdbuf()->pubsetbuf(buf,bufsize);
		ofs[i].rdbuf()->pubsetbuf(buf,bufsize);
	  }
	  if (!(ifs[0].good() && ifs[1].good() && ofs[0].good() && ofs[1].good())) {
		break;
	  }

	  logger.notice() << "Checking files: " << endl
		<< "  " << fn[0] << " " << fn[1];

	  size_t  pos, pos_end;
	  ifs[0].seekg(0, ios_base::end);
	  pos_end = ifs[0].tellg();
	  ifs[0].seekg(0, ios_base::beg);

	  vector<SEQ>  work_seq[2];

	  while (1) {  // read at _max_seq sequences.
		nlist n_list[2];
		work_seq[0].clear();
		work_seq[1].clear();
		if ((*f).first > FILETYPE::FASTA) {
		    const int qual_minus = (*f).first == FILETYPE::FASTQ33 ? FASTQ_PARAM::FASTQ_NEW : FASTQ_PARAM::FASTQ_OLD;
		    // read sequence and quality
		    read_fastq_ecc(ifs[0], work_seq[0], n_list[0], bind2nd(minus<char>(),(char)qual_minus));
		    read_fastq_ecc(ifs[1], work_seq[1], n_list[1], bind2nd(minus<char>(),(char)qual_minus));
		    is_fastq=true;
		} else {
		    read_fasta_ecc(ifs[0], work_seq[0], n_list[0]);
		    read_fasta_ecc(ifs[1], work_seq[1], n_list[1]);
		    is_fastq=false;
		}
		if (work_seq[0].empty() || work_seq[1].empty()) {
		    break;  // eof or error
		}
#       pragma omp parallel for schedule (dynamic) shared(n_list)
		for (i=0; i<work_seq[0].size(); i++) {
		    vector<char>  nk0, nk1;
		    trim_adapter(adapt, work_seq[0][i], work_seq[1][i],n_list[0][i],n_list[1][i], adapt_result[omp_get_thread_num()]);
		    const int trimmed = _trim_low_qual(work_seq[0][i],work_seq[1][i], n_list[0][i],n_list[1][i],qual_result[omp_get_thread_num()]) ;
		}
		_output(work_seq[0], work_seq[1], ofs[0], ofs[1], n_list[0], n_list[1], output_th_len,is_fastq);
		pos = (ifs[0].eof() ? pos_end : (size_t)ifs[0].tellg());
		logger.notice() << "\r  " << fn[0] << " " << fn[1] << "  (" << (int)(pos*100.0/pos_end) << "%)";
	  }
	  logger.notice() << "\r  " << fn[0] << " " << fn[1] << "  (100%)" << endl;
    }

    for (unsigned i = 1; i < NUM_OF_THREADS; ++i) {
	  adapt_result[0] += adapt_result[i];
	  qual_result[0] += qual_result[i];
    }

    logger.notice() << endl
	  << "Number of trimmed read with adapter: " << endl;
    logger.notice() << adapt_result[0] << endl;
    logger.notice() << endl
	  << "Number of trimmed read because of low quality or too short (< " << HIGH_QUALITY_REGION_TH << "bp): " << endl;
    logger.notice() << qual_result[0] << endl;

    return 0;
}


int Trim::trim_adapter(ADAPT_KMER& adapt, SEQ &seq_f, SEQ &seq_r,vector<int> &n_list1,vector<int> &n_list2, Result &res, bool is_fastq)
{
    int       i, j, k, shift;
    unsigned  fwd, rev, mask;
    unordered_map<unsigned,int>            fw;
    unordered_map<unsigned,int>::iterator  itr;
    vector<BLOCK> hit, ad[2];
    int offset = 0;
    BLOCK work, whit;
    if (seq_f.seq.size()<adapt.len || seq_r.seq.size()<adapt.len) {
	  return 0;
    }

    mask = (1<<adapt.len*2)-1;
    shift = (adapt.len-1)*2;
    ad[0].push_back (work);

    fwd = (unsigned)seq_f.seq[0];
    for (i = 1; i < adapt.len-1; i++) fwd = (fwd<<2)|(unsigned)seq_f.seq[i];
    for (; i < seq_f.seq.size (); i++) {
	  fwd = (fwd<<2)|(unsigned)seq_f.seq[i];
	  fwd &= mask;
	  fw[fwd] = i-adapt.len+1;

	  // find adaptor
	  itr = adapt.find(fwd);
	  if (itr == adapt.end()) {
		continue;
	  }
	  work.flag = 0;
	  if ((*itr).second > (adapt.OFFSET * 3) - (adapt.OFFSET >> 1))
		offset = adapt.OFFSET * 3;
	  else if ((*itr).second > (adapt.OFFSET * 2) - (adapt.OFFSET >> 1))
		offset = adapt.OFFSET * 2;
	  else if ((*itr).second > adapt.OFFSET - (adapt.OFFSET >> 1))
		offset = adapt.OFFSET;
	  else
		offset = 0;
	  k = (*itr).second - offset;
	  if (k > 0) {
		work.str = '+';
		work.sst = k+adapt.len-1;
		work.sed = k;
	  } else {
		work.str = '-';
		work.sst = k;
		work.sed = k-adapt.len+1;
	  }
	  work.qst = i-adapt.len+2;
	  work.qed = i+1;

	  BLOCK &b = ad[0].back();
	  if (b.flag != work.flag
		    || b.str != work.str
		    || b.sst <= work.sst
		    || abs (b-work) > 1) {
		ad[0].push_back (work);
	  } else {
		b.sed = work.sed;
		b.qed = work.qed;
	  }
    }

    work.flag = -1;
    ad[1].push_back (work); // dummy
    hit.push_back (work); // dummy

    rev = (unsigned)(seq_r.seq[0]^3)<<shift;
    for (i = 1; i < adapt.len-1; i++) rev = (rev>>2)|(unsigned)(seq_r.seq[i]^3)<<shift;
    for (; i < seq_r.seq.size (); i++) {
	  rev = (rev>>2)|(unsigned)(seq_r.seq[i]^3)<<shift;

	  itr = adapt.find(rev);
	  if (itr != adapt.end()) {
		work.flag = 0;
		if ((*itr).second > (adapt.OFFSET * 3) - (adapt.OFFSET >> 1))
		    offset = adapt.OFFSET * 3;
		else if ((*itr).second > (adapt.OFFSET * 2) - (adapt.OFFSET >> 1))
		    offset = adapt.OFFSET * 2;
		else if ((*itr).second > adapt.OFFSET - (adapt.OFFSET >> 1))
		    offset = adapt.OFFSET;
		else
		    offset = 0;
		k = (*itr).second - offset;
		if (k < 0) {
		    work.str = '+';
		    work.sst = -k+adapt.len-1;
		    work.sed = -k;
		} else {
		    work.str = '-';
		    work.sst = -k;
		    work.sed = -k-adapt.len+1;
		}
		work.qst = i-adapt.len+2;
		work.qed = i+1;

		BLOCK &b = ad[1].back();
		if (b.flag != work.flag
			  || b.str != work.str
			  || b.sst <= work.sst
			  || abs (b-work) > 1) {
		    ad[1].push_back (work);
		} else {
		    b.sed = work.sed;
		    b.qed = work.qed;
		}
	  }

	  // override if found in fowward matches.
	  itr = fw.find(rev);
	  if (itr != fw.end()) {
		work.sst = (*itr).second+adapt.len;
		work.qed = i+1;
		BLOCK &h = hit.back();
		if (h.sst<=work.sst || abs(h-work)>1) {
		    hit.push_back(work);
		} else {
		    h.sed = work.sed;
		    h.qed = work.qed;
		}
	  }
    }

    if (ad[0].size()==1 && ad[1].size()==1 && hit.size()==1) {
	  return 0;
    }
    if (hit.size() <= 1) {
	  hit.clear();
    } else {
	  for (i = 1; i < hit.size (); i++)
		for (j = i+2; j < hit.size (); j++)
		    if (hit[i].sst > hit[j].sst && abs(hit[i]-hit[j]) <= 1) {
			  hit[i].sed = hit[j].sed;
			  hit[i].qed = hit[j].qed;
			  hit[j].flag = 99;
		    }
	  k = 30;
	  for (i = 1; i < hit.size (); i++)
		if (hit[i].flag != 99 && hit[i].qed-hit[i].qst+1 > k) {
		    whit = hit[i];
		    k = hit[i].qed-hit[i].qst+1;
		}
	  hit.clear ();
	  if (k > 30) hit.push_back (whit);
    }

    int st[2] = {0, 0};
    int ed[2] = {static_cast<int>(seq_f.seq.size()), static_cast<int>(seq_r.seq.size())};
    char rep = 0;

    for (i=0; i<2; i++) {
	  SEQ &s = (i==0 ? seq_f : seq_r);
	  if (!is_fastq) {  // fasta
		st[i] = s.seq.size() - 1;
		ed[i] = -1;
	  // } else {
		// for (st[i]=0; st[i]<s.seq.size(); st[i]++) if (s.qual[st[i]]>2) break;
		// for (ed[i]=s.seq.size()-1; ed[i]>=0; ed[i]--) if (s.qual[ed[i]]>2) break;
	  }

	  if (ed[i]-st[i]+1 <= 30) {
		st[i] = ed[i] = seq_f.seq.size ();
	  } else {
		j = st[i]+1; k = 0;

		auto a=ad[i].begin(), a_e=ad[i].end();
		for (a++; a!=a_e; a++) {
		    if ((*a).qst > j) j = (*a).qst;
		    if ((*a).qed > ed[i]+1) { k += ed[i]+1-j+1; break; }
		    else k += (*a).qed-j+1;
		    j = (*a).qed+1;
		}

		if (ed[i]-st[i]+1-k < 30) {
		    rep = 1;
		} else {
		    for (a=ad[i].begin()+1,j=0; a!=a_e; a++) {
			  if ((*a).qst-st[i] >= 16) continue;
			  if ((*a).str == '+' && (*a).sed < 16) j = (*a).qed+(*a).sed-1;
			  else j = (*a).qed;
		    }
		    if (j > st[i]) st[i] = j;
		}
	  }
    }

    if (rep) {
	  st[0] = ed[0] = seq_f.seq.size();
	  st[1] = ed[1] = seq_r.seq.size();
    } else {
	  vector<int>::iterator it = n_list1.begin(),it_ed=n_list1.end();
	  for (i=0; i<2; i++) {
		auto a=ad[i].begin(), a_e=ad[i].end();
		for (a++,j=9999; a!=a_e; a++) {
		    if ((*a).str == '-') {
			  if ((*a).qed > ed[i]-15 || 
				    ((*a).sst < 5 && (*a).qed-(*a).qst+1 >= 17)) {
				j = (*a).qst-(*a).sst-1;
				break;
			  }
		    }
		}
		if (j < ed[i]) ed[i] = j;
		for (; it != it_ed; ++it) {
		    *it -= st[i];
		}
		it = n_list2.begin();
		it_ed = n_list2.end();
	  }
	  if (!hit.empty()) {
		BLOCK &h = hit[0];
		int w1 = h.sst+h.qst-1;
		int w2 = h.qed+h.sed-1;
		if (w1 == w2
			  && w1 < seq_f.seq.size ()
			  && h.sed < 16
			  && h.qst < 16) {
		    if (ad[0].size () == 1
				&& ad[1].size () == 1
				&& (ed[0]+1-h.sst >= 40 || ed[1]+1-h.qed >= 40)) {
			  ;
		    } else {
			  if (h.sst-1 < ed[0]) ed[0] = h.sst-1;
			  if (h.qed-1 < ed[1]) ed[1] = h.qed-1;
		    }
		}
	  }
    }

    // set 'N' at adaptor (4 is 'N')
    // vector<char>::iterator s, s_e, q;
    // if (!is_fastq) {  // FASTA
    // for (s=seq_f.seq.begin(),s_e=s+st[0]; s!=s_e; s++) { *s=4; ++res.num_base_f; }
    // for (s+=(ed[0]-st[0]+1),s_e=seq_f.seq.end(); s<s_e; s++) { *s=4; ++res.num_base_f; }
    // for (s=seq_r.seq.begin(),s_e=s+st[1]; s!=s_e; s++) { *s=4; ++res.num_base_r; }
    // for (s+=(ed[1]-st[1]+1),s_e=seq_r.seq.end(); s<s_e; s++) { *s=4; ++res.num_base_r; }
    // } else {  // FASTQ
    // unsigned long long *base_ptr;
    // #define ASSIGN_SEQ(FR,N)                               \
    // s = seq_##FR.seq.begin();                          \
    // q = seq_##FR.qual.begin();                         \
    // base_ptr = &(res.num_base_##FR);                     \
    // s_e = std::min(s+st[N], seq_##FR.seq.end());       \
    // for (; s!=s_e; s++,q++) { *s=4; *q=Q_OPTION_LOWER_TH; ++(*base_ptr); }            \
    // s_e = seq_##FR.seq.end();                          \
    // if (s!=s_e) {                                      \
    // i = std::max(std::min(ed[N]-st[N]+1, (int)(s_e-s)), 0);       \
    // for (s+=i,q+=i; s<s_e; s++,q++) { *s=4; *q=Q_OPTION_LOWER_TH; ++(*base_ptr); }  \
    // }
    // ASSIGN_SEQ(f,0);
    // ASSIGN_SEQ(r,1);
    // #undef ASSIGN_SEQ
    // }
    // if ((st[0]>0 || ed[0]<seq_f.seq.size() - 1) && (st[1]>0 || ed[1]<seq_r.seq.size() - 1)) {
    // ++res.num_read_r;
    // ++res.num_read_f;
    // ++res.num_pair_and;
    // ++res.num_pair_or;
    // } else if (st[0]>0 || ed[0]<seq_f.seq.size() - 1) {
    // ++res.num_read_f;
    // ++res.num_pair_or;
    // } else if (st[1]>0 || ed[1]<seq_r.seq.size() - 1) {
    // ++res.num_read_r;
    // ++res.num_pair_or;
    // }
    // return st[0]>0 || st[1]>0 || ed[0]<seq_f.seq.size() - 1 || ed[1]<seq_r.seq.size() - 1 ? 1 : 0;
    vector<char>::iterator s, s_e, q;
    unsigned trimmed_num = 0;
    if (!is_fastq) {  // FASTA
	  SEQ new_f, new_r;
	  if (st[0] != 0 || ed[0] != seq_f.seq.size()-1) {
		for (s=seq_f.seq.begin()+st[0], s_e=std::min(seq_f.seq.begin()+ed[0]+1, seq_f.seq.end()); s<s_e; ++s) {
		    new_f.seq.push_back(*s);
		}
		res.num_base_f += seq_f.seq.size() - new_f.seq.size();
		seq_f.seq.swap(new_f.seq);
		++res.num_read_f;
		++trimmed_num;
	  }
	  if (st[1] != 0 || ed[1] != seq_r.seq.size()-1) {
		for (s=seq_r.seq.begin()+st[1], s_e=std::min(seq_r.seq.begin()+ed[1]+1, seq_r.seq.end()); s<s_e; ++s) {
		    new_r.seq.push_back(*s);
		}
		res.num_base_r += seq_r.seq.size() - new_r.seq.size();
		seq_r.seq.swap(new_r.seq);
		++res.num_read_r;
		++trimmed_num;
	  }
    } else {  // FASTQ
	  SEQ new_f, new_r;
	  if (st[0] != 0 || ed[0] != seq_f.seq.size()-1) {
		for (s=seq_f.seq.begin()+st[0], q=seq_f.qual.begin()+st[0], s_e=std::min(seq_f.seq.begin()+ed[0]+1, seq_f.seq.end()); s<s_e; ++s, ++q) {
		    new_f.seq.push_back(*s);
		    new_f.qual.push_back(*q);
		}
		res.num_base_f += seq_f.seq.size() - new_f.seq.size();
		seq_f.seq.swap(new_f.seq);
		seq_f.qual.swap(new_f.qual);
		++res.num_read_f;
		++trimmed_num;
	  }
	  if (st[1] != 0 || ed[1] != seq_r.seq.size()-1) {
		for (s=seq_r.seq.begin()+st[1], q=seq_r.qual.begin()+st[1], s_e=std::min(seq_r.seq.begin()+ed[1]+1, seq_r.seq.end()); s<s_e; ++s, ++q) {
		    new_r.seq.push_back(*s);
		    new_r.qual.push_back(*q);
		}
		res.num_base_r += seq_r.seq.size() - new_r.seq.size();
		seq_r.seq.swap(new_r.seq);
		seq_r.qual.swap(new_r.qual);
		++res.num_read_r;
		++trimmed_num;
	  }
    }

    if (trimmed_num == 2) {
	  ++res.num_pair_and;
    }

    if (trimmed_num >= 1) {
	  ++res.num_pair_or;
    }

    return trimmed_num > 0 ? 1 : 0;
}


int Trim::_exec_internal(ifstream& kfs, Trim::ADAPT_KMER &adaptn, Trim::INTERNAL_ADAPT& adapt)
{

    omp_set_num_threads(NUM_OF_THREADS);

    // load k-mers and put partial on memory.

    vector<Result> internal_result(NUM_OF_THREADS);
    vector<Result> adapt_result(NUM_OF_THREADS);
    vector<Result> qual_result(NUM_OF_THREADS);
    vector<Result> tmp_result(NUM_OF_THREADS);
    INTERNAL_ADAPT first_half;
    INTERNAL_ADAPT second_half;
    bool is_fastq;
    int output_th_len=atoi(exec_option['l'].c_str());

    if(atoi(exec_option['a'].c_str())==2){
	  first_half=INTERNAL_ADAPT ("CTGTCTCTTATACACATCT",0);
	  second_half=INTERNAL_ADAPT ("AGATGTGTATAAGAGACAG",0);
    }

    for (auto f=input_file.begin(),f_e=input_file.end(); f!=f_e; f+=2) {
	  // generate output file pathes
	  ofstream dbg;
	  string fn[4]={(*f).second, (*(f+1)).second};
	  fn[2] = fn[0] + ".int_trimmed";
	  fn[3] = fn[1] + ".int_trimmed";

	  // open files
	  int i;
	  ifstream ifs[2];
	  ofstream ofs[2];
	  const int bufsize=1000000;
	  char buf[bufsize];
	  for (i=0; i<2; i++) {
		ifs[i].open(fn[i].c_str());
		ofs[i].open(fn[i+2].c_str());
		if (!ifs[i].good()) {
		    logger.error() << "ERROR: cannot open: " << fn[i] << endl;
		}
		if (!ofs[i].good()) {
		    logger.error() << "ERROR: cannot open: " << fn[i+2] << endl;
		}
		ifs[i].rdbuf()->pubsetbuf(buf,bufsize);
		ofs[i].rdbuf()->pubsetbuf(buf,bufsize);
	  }
	  if (!(ifs[0].good() && ifs[1].good() && ofs[0].good() && ofs[1].good())) {
		break;
	  }

	  logger.notice() << "Checking files: " << endl
		<< "  " << fn[0] << " " << fn[1];

	  size_t  pos, pos_end;
	  ifs[0].seekg(0, ios_base::end);
	  pos_end = ifs[0].tellg();
	  ifs[0].seekg(0, ios_base::beg);

	  vector<SEQ>  work_seq[2];

	  while (1) {  // read at _max_seq sequences.
		nlist n_list[2];
		work_seq[0].clear();
		work_seq[1].clear();
		if ((*f).first > FILETYPE::FASTA) {
		    const int qual_minus = (*f).first == FILETYPE::FASTQ33 ? FASTQ_PARAM::FASTQ_NEW : FASTQ_PARAM::FASTQ_OLD;
		    // read sequence and quality
		    read_fastq_ecc(ifs[0], work_seq[0], n_list[0], bind2nd(minus<char>(),(char)qual_minus));
		    read_fastq_ecc(ifs[1], work_seq[1], n_list[1], bind2nd(minus<char>(),(char)qual_minus));
		    is_fastq=true;
		}else{
		    read_fasta_ecc(ifs[0], work_seq[0], n_list[0]);
		    read_fasta_ecc(ifs[1], work_seq[1], n_list[1]);
		    is_fastq=false;
		}
		if (work_seq[0].empty() || work_seq[1].empty()) {
		    break;  // eof or error
		}
#       pragma omp parallel for schedule (dynamic)  shared(n_list)
		for (i=0; i<work_seq[0].size(); i++) {
		    int is_trimmed = trim_internal(adapt, work_seq[0][i], work_seq[1][i], internal_result[omp_get_thread_num()],is_fastq);
		    if (!is_trimmed) {
			  is_trimmed+=trim_internal(first_half, work_seq[0][i], work_seq[1][i], tmp_result[omp_get_thread_num()],is_fastq);
			  is_trimmed+=trim_internal(second_half, work_seq[0][i], work_seq[1][i], tmp_result[omp_get_thread_num()],is_fastq);
			  if (!is_trimmed && exec_option['r'] == "false") {
				trim_adapter(adaptn, work_seq[0][i], work_seq[1][i],n_list[0][i], n_list[1][i], adapt_result[omp_get_thread_num()],is_fastq);
			  }
		    }

		    if ((*f).first > FILETYPE::FASTA) {
			  is_trimmed = _trim_low_qual(work_seq[0][i],work_seq[1][i], n_list[0][i],n_list[1][i],qual_result[omp_get_thread_num()]) ;
		    }
		}
		_output(work_seq[0], work_seq[1], ofs[0], ofs[1], n_list[0], n_list[1], output_th_len,is_fastq);
		if ((*f).first > FILETYPE::FASTA)
		    pos = (ifs[0].eof() ? pos_end : (size_t)ifs[0].tellg());
		logger.notice() << "\r  " << fn[0] << " " << fn[1] << "  (" << (int)(pos*100.0/pos_end) << "%)";
	  }
	  logger.notice() << "\r  " << fn[0] << " " << fn[1] << "  (100%)" << endl;
    }

    for (unsigned i = 1; i < NUM_OF_THREADS; ++i) {
	  internal_result[0] += internal_result[i];
	  adapt_result[0] += adapt_result[i];
	  qual_result[0] += qual_result[i];
	  tmp_result[0] += tmp_result[i];
    }
    logger.notice() << endl
	  << "Number of trimmed read with internal adapter: " << endl;
    logger.notice() << internal_result[0] << endl;
    logger.notice() << endl
	  << "Number of trimmed read with adapter: " << endl;
    logger.notice() << adapt_result[0] << endl;
    logger.notice() << endl
	  << "Number of trimmed read because of low quality or too short (< " << HIGH_QUALITY_REGION_TH << "bp): " << endl;
    logger.notice() << qual_result[0] << endl;
    // logger.notice() << endl
	  // << "half forward" << endl;
    // logger.notice() << tmp_result[0] << endl;

    return 0;
}


int Trim::trim_internal(const Trim::INTERNAL_ADAPT& adapt, SEQ &seq_f, SEQ &seq_r, Result &res, bool is_fastq)
{
    int i = 0;
    int code = 0;
    vector<char> fwd, rev;
    int st[2] = {0, 0};
    int ed[2] = {static_cast<int>(seq_f.seq.size()), static_cast<int>(seq_r.seq.size())};
    int for_end;

    if (seq_f.seq.size() >= adapt.len) {
	  for_end = min(adapt.all_len, static_cast<int>(seq_f.seq.size()));
	  fwd.resize(adapt.len);
	  for (i = 0; i < adapt.len; ++i) {
		fwd[i] = seq_f.seq[i];
	  }

	  auto it = adapt.find(fwd);
	  if (it != adapt.end()) {
		st[0] = seq_f.seq.size() - 1;
		code = it->second;
	  } else {
		for (i = adapt.len - adapt.all_len +1; i < (int)seq_f.seq.size() - adapt.len + 1; ++i) {
		    int buf=i-1;
		    if (i > 0){
			  fwd.erase(fwd.begin());
		    }
		    else buf=0;
		    if (i + adapt.all_len - 1 < (int)seq_f.seq.size()) {
			  fwd.push_back(seq_f.seq[i + adapt.all_len - 1]);
		    }
		    auto it = adapt.find(fwd);
		    if (it != adapt.end()) {
			  ed[0] = buf;
			  code = it->second;
			  break;
		    }
		}
	  }
    }


    if (seq_r.seq.size() >= adapt.len) {
	  for_end = min(adapt.all_len, static_cast<int>(seq_r.seq.size()));
	  rev.resize(adapt.len);
	  for (i = 0; i < adapt.len; ++i) {
		rev[i] = seq_r.seq[i];
	  }
	  auto it = adapt.find(rev);
	  if (it != adapt.end()) {
		//if (code * it->second <= 0) {
		st[1] = seq_r.seq.size() - 1;
		code = it->second;
		//} else {
		//  ed[0] = 0;
		//  st[1] = seq_r.seq.size() - 1;
		//}
	  } else {
		for (i = adapt.len - adapt.all_len +1 ; i < (int)seq_r.seq.size() - adapt.len + 1; ++i) {
		    int buf=i-1;
		    if (i > 0){
			  rev.erase(rev.begin());
		    }
		    else buf=0;
		    if (i + adapt.all_len - 1 < (int)seq_r.seq.size()) {
			  rev.push_back(seq_r.seq[i + adapt.all_len - 1]);
		    }
		    auto it = adapt.find(rev);
		    if (it != adapt.end()) {
			  //if (code * it->second <= 0) {
			  ed[1] = buf;
			  code = it->second;
			  //} else {
			  //  ed[0] = 0;
			  //  st[1] = seq_r.seq.size() - 1;
			  //}
			  break;
		    }
		}
	  }
    }

    vector<char>::iterator s, s_e, q;
    unsigned trimmed_num = 0;
    if (!is_fastq) {  // FASTA
	  SEQ new_f, new_r;
	  if (st[0] != 0 || ed[0] != seq_f.seq.size()) {
		for (s=seq_f.seq.begin()+st[0], s_e=std::min(seq_f.seq.begin()+ed[0]+1, seq_f.seq.end()); s<s_e; ++s) {
		    new_f.seq.push_back(*s);
		}
		res.num_base_f += seq_f.seq.size() - new_f.seq.size();
		seq_f.seq.swap(new_f.seq);
		++res.num_read_f;
		++trimmed_num;
	  }
	  if (st[1] != 0 || ed[1] != seq_r.seq.size()) {
		for (s=seq_r.seq.begin()+st[1], s_e=std::min(seq_r.seq.begin()+ed[1]+1, seq_r.seq.end()); s<s_e; ++s) {
		    new_r.seq.push_back(*s);
		}
		res.num_base_r += seq_r.seq.size() - new_r.seq.size();
		seq_r.seq.swap(new_r.seq);
		++res.num_read_r;
		++trimmed_num;
	  }
    } else {  // FASTQ
	  SEQ new_f, new_r;
	  if (st[0] != 0 || ed[0] != seq_f.seq.size()) {
		for (s=seq_f.seq.begin()+st[0], q=seq_f.qual.begin()+st[0], s_e=std::min(seq_f.seq.begin()+ed[0]+1, seq_f.seq.end()); s<s_e; ++s, ++q) {
		    new_f.seq.push_back(*s);
		    new_f.qual.push_back(*q);
		}
		res.num_base_f += seq_f.seq.size() - new_f.seq.size();
		seq_f.seq.swap(new_f.seq);
		seq_f.qual.swap(new_f.qual);
		++res.num_read_f;
		++trimmed_num;
	  }
	  if (st[1] != 0 || ed[1] != seq_r.seq.size()) {
		for (s=seq_r.seq.begin()+st[1], q=seq_r.qual.begin()+st[1], s_e=std::min(seq_r.seq.begin()+ed[1]+1, seq_r.seq.end()); s<s_e; ++s, ++q) {
		    new_r.seq.push_back(*s);
		    new_r.qual.push_back(*q);
		}
		res.num_base_r += seq_r.seq.size() - new_r.seq.size();
		seq_r.seq.swap(new_r.seq);
		seq_r.qual.swap(new_r.qual);
		++res.num_read_r;
		++trimmed_num;
	  }
    }

    if (trimmed_num == 2) {
	  ++res.num_pair_and;
    }

    if (trimmed_num >= 1) {
	  ++res.num_pair_or;
    }

    /*
#define ASSIGN_SEQ(FR,N)                               \
s = seq_##FR.seq.begin();                          \
q = seq_##FR.qual.begin();                         \
s_e = std::min(s+st[N], seq_##FR.seq.end());       \
for (; s!=s_e; s++,q++) { *s=4; *q=0; }            \
s_e = seq_##FR.seq.end();                          \
if (s!=s_e) {                                      \
i = std::min(ed[N]-st[N]+1, (int)(s_e-s));       \
for (s+=i,q+=i; s<s_e; s++,q++) { *s=4; *q=0; }  \
}
ASSIGN_SEQ(f,0);
ASSIGN_SEQ(r,1);
#undef ASSIGN_SEQ
}
SEQ new_f;
for (auto its = seq_f.seq.begin(), itq = seq_f.qual.begin(), ends = seq_f.seq.end(); its != ends; ++its, ++itq) {
if (*its == 4 && *itq == 0) { break; }
new_f.seq.push_back(*its);
new_f.qual.push_back(*itq);
}
seq_f.seq.swap(new_f.seq);
seq_f.qual.swap(new_f.qual);

SEQ new_r;
for (auto its = seq_r.seq.begin(), itq = seq_r.qual.begin(), ends = seq_r.seq.end(); its != ends; ++its, ++itq) {
if (*its == 4 && *itq == 0) { break; }
new_r.seq.push_back(*its);
new_r.qual.push_back(*itq);
}
seq_r.seq.swap(new_r.seq);
seq_r.qual.swap(new_r.qual);
*/
return trimmed_num > 0 ? 1 : 0;
}

int Trim::_trim_low_qual(SEQ &seq_f,SEQ &seq_r, vector<int> &n_list1,vector<int> &n_list2,Result &res)
{
    SEQ *seq=&seq_f;
    int trimmed_f=0;
    int trimmed_r=0;
    int *trimmed=&trimmed_f;
    vector<int> *n_list=&n_list1;

    for(int k=2;k>0;--k){
	  int i, j, st, ed, l;
	  SEQ work;

	  st = 0;
	  ed = seq->seq.size ();

	  for (i=0,l=0; i<ed; i++) {
		if (seq->qual[i] > th_qv) {
		    l++;
		    if (l == HIGH_QUALITY_REGION_TH) break;
		}else {
		    st += l+1;
		    l = 0;
		}
	  }
	  for (i=ed-1,l=0; i>=0; i--) {
		if (seq->qual[i] > th_qv) {
		    l++;
		    if (l == HIGH_QUALITY_REGION_TH) break;
		}else {
		    ed -= l+1;
		    l = 0;
		}
	  }

	  vector<int>::iterator it = n_list->begin();
	  for (; it != n_list->end(); ++it) {
		*it -= st;
	  }

	  // if (!nk.empty()) {
	  // l = seq.seq.size() - nk.size();
	  // for (i = st; i < ed-l; i++) {
	  // if (nk[i] == 0) {
	  // unsigned size = seq.seq.size();
	  // seq.seq.clear();
	  // trimmed= size;
	  // }
	  // }
	  // }

	  for (i = st; i < ed; i++) work.seq.push_back (seq->seq[i]);
	  seq->seq.swap (work.seq);
	  if (!seq->qual.empty()) {
		for (i = st; i < ed; i++) work.qual.push_back (seq->qual[i]);
		seq->qual.swap (work.qual);
	  }
	  if (st < ed)
		*trimmed =st + work.seq.size() - ed;
	  else
		*trimmed= work.seq.size();

	  seq=&seq_r;
	  trimmed=&trimmed_r;
	  n_list=&n_list2;
    }
    res.num_base_f += trimmed_f;
    res.num_read_f += trimmed_f > 0 ? 1 : 0;
    res.num_base_r += trimmed_r;
    res.num_read_r += trimmed_r > 0 ? 1 : 0;
    // res.num_pair_and += (trimmed_f > 0 ? 1 : 0 + trimmed_r > 0 ? 1 : 0) >> 1;
    res.num_pair_and += (trimmed_f * trimmed_r) > 0 ? 1 : 0;
    res.num_pair_or += (trimmed_f + trimmed_r) > 0 ? 1 : 0;

    return res.num_read_f*res.num_read_r > 0 ? 1 : 0;
}




void Trim::_output(vector<SEQ> &seq1, vector<SEQ> &seq2, ofstream &ofsp, ofstream &ofss, nlist &n_list1, nlist &n_list2, int len,bool is_fastq)
{
    auto s1=seq1.begin(), s1_e=seq1.end(), s2=seq2.begin();
    if (len>=0) {
	  for (int i = 0; s1!=s1_e; s1++,s2++,i++) {
		if ((*s1).seq.size () >= len) {
		    if ((*s2).seq.size () >= len) {
			  (*s1).put(ofsp, n_list1[i],is_fastq);
			  (*s2).put(ofss, n_list2[i],is_fastq);
		    } else {
			  (*s2).seq.clear();
			  (*s2).qual.clear();
			  if(exec_option['f']=="false"){
				(*s1).seq.clear();
				(*s1).qual.clear();
				continue;
			  }
			  (*s1).put(ofsp, n_list1[i],is_fastq);
			  (*s2).put(ofss, n_list2[i],is_fastq);
		    }
		} else if ((*s2).seq.size () >= len) {
		    (*s1).seq.clear();
		    (*s1).qual.clear();
		    if(exec_option['f']=="false"){
			  (*s2).seq.clear();
			  (*s2).qual.clear();
			  continue;
		    }
		    (*s1).put(ofsp, n_list1[i],is_fastq);
		    (*s2).put(ofss, n_list2[i],is_fastq);
		}
	  }
    } else {
	  // only trim adapter mode
	  for (int i = 0; s1!=s1_e; s1++,s2++,i++) {
		(*s1).put(ofsp, n_list1[i], is_fastq);
		(*s2).put(ofss, n_list2[i], is_fastq);
	  }
    }
}

