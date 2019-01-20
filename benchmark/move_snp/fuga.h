#ifndef FUGA_H_
#define FUGA_H_

#include "hoge.h"

class Fuga;
class Aln;
class RefPosInfo;
class QryPosInfo;

class Aln {
 public:
  Aln(const int64_t ref_beg, const int64_t qry_beg, const int64_t ref_end,
      const int64_t qry_end, const std::string& qry_id, const char qry_strand)
      : ref_beg_(ref_beg),
        qry_beg_(qry_beg),
        ref_len_(ref_end - ref_beg + 1),
        qry_len_(std::abs(qry_end - qry_beg) + 1),
        qry_id_(qry_id),
        qry_strand_(qry_strand),
        ref_aln_seq_(""),
        qry_aln_seq_("") {}
  ~Aln() {}
  const std::string& qry_id() const { return qry_id_; }
  const char qry_strand() const { return qry_strand_; }
  int64_t RefPos2QryPos(int64_t ref_pos) const;
  void SetSeq(int64_t i);
 
 private:
  int64_t ref_beg_;
  int64_t qry_beg_;
  int64_t ref_len_;
  int64_t qry_len_;
  std::string qry_id_;
  char qry_strand_;
  std::string ref_aln_seq_;
  std::string qry_aln_seq_;
  int64_t RefPos2AlnPos(const int64_t ref_pos) const;
  int64_t AlnPos2QryPos(const int64_t aln_pos) const;
};

class RefPosInfo {
 public:
  RefPosInfo() : depth_(0), aln_idx_(0) {}
  ~RefPosInfo() {}
  int depth() const { return depth_; }
  int aln_idx() const { return aln_idx_; }
  void SetInfo(const int64_t aln_idx) {
    aln_idx_ = aln_idx;
    ++depth_;
  }

 private:
  int depth_;
  int aln_idx_;
};

class QryPosInfo {
 public:
  QryPosInfo() : depth_(0) {}
  ~QryPosInfo() {}
  int depth() const { return depth_; }
  void SetInfo() { ++depth_; }

 private:
  int depth_;
};

class Fuga {
 public:
  Fuga(string RootAlnFNom, string IsoAlnFNom, string DeltaFNom, string RefGenomeFNom, string RootGenomeFNom, string OutDir) :
    kRootAlnFNom_(RootAlnFNom), kIsoAlnFNom_(IsoAlnFNom),
    kDeltaFNom_(DeltaFNom), kRefGenomeFNom_(RefGenomeFNom), kRootGenomeFNom_(RootGenomeFNom), kOutDir_(OutDir),
    kRandSeed_(1), kUniqAlnLen_(1000), AlnLen_(0), IsoNum_(0), SnpNum_(0), root_size_(0) {}
  ~Fuga() {}
  void Main();
  u32 Pos2AlnPos(u32 Pos) const { return Pos2AlnPos_.at(Pos); }
  u32 AlnPos2Pos(u32 AlnPos) const { return AlnPos2Pos_.at(AlnPos); }

 private:
  const string kRootAlnFNom_;
  const string kIsoAlnFNom_;
  const string kDeltaFNom_;
  const string kRefGenomeFNom_;
  const string kRootGenomeFNom_;
  const string kOutDir_;
  const u32 kRandSeed_;
  const u32 kUniqAlnLen_;
  void SetAln();
  void SetOrgSnp();
  void ReadDelta();
  void MoveSnp();
  void ReadRoot();
  string RootSeq_;
  vector<u32> Pos2AlnPos_;
  vector<u32> AlnPos2Pos_;
  u32 AlnLen_;
  vector<string> Iso_;
  vector<string> IsoSeq_;
  u32 IsoNum_;
  u32 SnpNum_;
  vector<u32> OrgSnpPos_;
  vector<u32> AltSnpPos_;
  vector<vector<bool> > RootOrAlt_;

  std::vector<std::string> root_cntg_;
  std::unordered_map<std::string, std::string> root_;
  std::string root_merged_seq_;
  u32 root_size_;
  void add_2_aln(const int64_t ref_beg, const int64_t qry_beg,
                 const int64_t ref_end, const int64_t qry_end,
                 const std::string& qry_id, const char qry_strand) {
    aln_.push_back(Aln(ref_beg, qry_beg, ref_end, qry_end, qry_id, qry_strand));
  }
  bool JdgUniq(const std::string& contig_id, const int64_t pos, const char r_or_q) const;
  pair<string, u32> GetUniqPos();
  char GetAltAllele(char root_allele);
  char GetRevComp(char nuc) const;
  std::vector<Aln> aln_;
  std::map<std::string, std::vector<RefPosInfo>> ref_id_2_pos_info_;
  std::map<std::string, std::vector<QryPosInfo>> qry_id_2_pos_info_;
};

#endif  // FUGA_H_
