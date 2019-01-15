// -*- mode:c++ -*-
#ifndef __BINSTR_H__
#define __BINSTR_H__

////////////////////////////////////////////////////////////////////////////////
// 64bit

typedef unsigned long long  u64_t;
const u64_t U64MINUS1 = 0xFFFFFFFFFFFFFFFFull;


////////////////////////////////////////////////////////////////////////////////
// binary string (1char=2bit)

struct binstr_t {
  u64_t               *value;
  const unsigned short len;  // string length

  //////////////////////////////////////////////////////////////////////////////
  // constructors and destructor

  binstr_t():value(NULL),len(0) {}
  binstr_t(unsigned short l):len(l) { value=new u64_t[(len+31)/32]; clear(); }
  binstr_t(const binstr_t& rhs):value(NULL),len(rhs.len) { (*this)=rhs; }
  virtual ~binstr_t() { if(value) { delete[] value; } }


  //////////////////////////////////////////////////////////////////////////////
  // utilities

  inline void resize(unsigned short _len) {
    if (value) {
      delete[] value;
    }
    *(const_cast<unsigned short*>(&(this->len))) = _len;
    value = (_len>0 ? new u64_t[(_len+31)/32] : NULL);
    clear();
  }

  inline u64_t to_u64() const {
    u64_t *v=value;
    for (int i=1,n=(len+31)/32; i<n; i++) {
      if (*(++v)) {
        return U64MINUS1;
      }
    }
    return *value;
  }

  inline binstr_t& flip() {
    u64_t *v=value;
    for (int i=0,n=(len+31)/32; i<n; i++,v++) {
      *v ^= U64MINUS1;
    }
    if (len%32>0) {
      (*(--v)) &= ((0x1ull<<(2*(len%32)))-0x1ull);  // mask
    }
    return *this;
  }

  // character base (not bit)
  inline void set(int i, unsigned char c) {
    u64_t *v = value;
    int   j = i*2;
    while (j>=64) {
      j-=64;
      v++;
    }
    ((*v)&=(U64MINUS1^(0x3ull<<j))) |= ((c&0x3ull)<<j);
  }

  inline void clear() {
    u64_t *v=value;
    for (int i=0,n=(len+31)/32; i<n; i++) {
      *v++ = 0;
    }
  }

  //////////////////////////////////////////////////////////////////////////////
  // operators

  template<typename T>
  inline binstr_t& operator=(const T rhs) {
    if (!value) {
      resize(1);
    }
    u64_t *v=value;
    for (int i=1,n=(len+31)/32; i<n; i++) {
      *(++v) = 0;
    }
    *value = (len>31 ? rhs : (rhs & ((0x1ull<<(2*(len%32)))-0x1ull)));
    return *this;
  }

  inline binstr_t& operator=(const binstr_t& rhs) {
    if (len!=rhs.len || !value) {
      resize(rhs.len);
    }
    u64_t *v=value, *vr=rhs.value;
    for (int i=0,n=(len+31)/32; i<n; i++) {
      *v++ = *vr++;
    }
    return *this;
  }

  template<typename T>
  inline binstr_t operator&(T rhs) const {
    binstr_t b(len);
    u64_t *v=b.value;
    *v++ = (*value) & rhs;
    for (int i=1,n=(len+31)/32; i<n; i++) {
      *v++ = 0;
    }
    return b;
  }

  template<typename T>
  inline binstr_t operator|(T rhs) const {
    binstr_t b(*this); if (value) { *(b.value)|=rhs; } return b;
  }

  template<typename T>
  inline binstr_t operator^(T rhs) const {
    binstr_t b(*this); if (value) { *(b.value)^=rhs; } return b;
  }

  template<typename T>
  inline binstr_t& operator&=(T rhs) {
    u64_t *v=value;
    *v++ &= rhs;
    for (int i=1,n=(len+31)/32; i<n; i++) {
      *v++ = 0;
    }
    return *this;
  }    

  template<typename T>
  inline binstr_t& operator|=(T rhs) {
    if (value) { *value|=rhs; } return *this;
  }

  template<typename T>
  inline binstr_t& operator^=(T rhs) {
    if (value) { *value^=rhs; } return *this;
  }

#define DEFINE_OPERATOR(B,OPE)                                 \
    u64_t *v=(B).value, *vr=rhs.value;                         \
    int i, n=(len+31)/32, nr=(rhs.len+31)/32;                  \
    for (i=0; i<n && i<nr; i++) {                              \
      *v++ OPE *vr++;                                          \
    }                                                          \
    if (i==n && len%32>0) {                                    \
      (*(--v)) &= ((0x1ull<<(2*(len%32)))-0x1ull);  /* mask */ \
    }
  inline binstr_t operator&(const binstr_t& rhs) const {
    binstr_t b(*this); DEFINE_OPERATOR(b,&=); return b;
  }
  inline binstr_t operator|(const binstr_t& rhs) const {
    binstr_t b(*this); DEFINE_OPERATOR(b,|=); return b;
  }
  inline binstr_t operator^(const binstr_t& rhs) const {
    binstr_t b(*this); DEFINE_OPERATOR(b,^=); return b;
  }
  inline binstr_t& operator&=(const binstr_t& rhs) {
    DEFINE_OPERATOR(*this, &=); return *this;
  }
  inline binstr_t& operator|=(const binstr_t& rhs) {
    DEFINE_OPERATOR(*this, |=); return *this;
  }
  inline binstr_t& operator^=(const binstr_t& rhs) {
    DEFINE_OPERATOR(*this, ^=); return *this;
  }
#undef DEFINE_OPERATOR

#define DEFINE_OPERATOR(B)                                              \
    if (len*2<=n || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      int i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value;                                               \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v++ = *(v+nd);                                               \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v++ = ((*(v+nd))>>nm) | ((*(v+nd+1))<<(64-nm));              \
        }                                                               \
        *v++ = (*(v+nd))>>nm;                                           \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v++ = 0;                                                       \
      }                                                                 \
    }
  inline binstr_t operator>>(int n) const {
    binstr_t b(*this); DEFINE_OPERATOR(b); return b;
  }
  inline binstr_t& operator>>=(int n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

#define DEFINE_OPERATOR(B)                                              \
    if (len*2<=n || len==0) {                                           \
      (B).clear();                                                      \
    } else {                                                            \
      int i, nb=(len+31)/32, nd=n/64, nm=n%64;                          \
      u64_t *v=(B).value + nb - 1;                                      \
      if (nm==0) {                                                      \
        for (i=0; i<nb-nd; i++) {                                       \
          *v-- = *(v-nd);                                               \
        }                                                               \
      } else {                                                          \
        for (i=0; i<nb-nd-1; i++) {                                     \
          *v-- = ((*(v-nd))<<nm) | ((*(v-nd-1))>>(64-nm));              \
        }                                                               \
        *v-- = (*(v-nd))<<nm;                                           \
        i++;                                                            \
      }                                                                 \
      for (; i<nb; i++) {                                               \
        *v-- = 0;                                                       \
      }                                                                 \
      if (len%32>0) {                                                   \
        *((B).value+nb-1) &= ((0x1ull<<(2*(len%32)))-0x1ull);  /* mask */ \
      }                                                                 \
    }
  inline binstr_t operator<<(int n) const {
    binstr_t b(*this); DEFINE_OPERATOR(b); return b;
  }
  inline binstr_t& operator<<=(int n) {
    DEFINE_OPERATOR(*this); return *this;
  }
#undef DEFINE_OPERATOR

  inline bool operator==(const binstr_t& rhs) const {
    u64_t *v=value, *vr=rhs.value;
    int i, n=(len+31)/32, nr=(rhs.len+31)/32;
    for (i=0; i<n && i<nr; i++) {
      if (*v++ != *vr++) {
        return false;
      }
    }
    if (n<nr) {
      for (; i<nr; i++) {
        if (*vr++) {
          return false;
        }
      }
    } else {
      for (; i<n; i++) {
        if (*v++) {
          return false;
        }
      }
    }
    return true;
  }

  inline bool operator<(const binstr_t& rhs) const {
    int i=(len+31)/32-1;
    u64_t *l=value, *r=rhs.value;
    while (i && *(l+i)==*(r+i)) { i--; }
    return (*(l+i)<*(r+i));
  };

  // character base (not bit)
  inline unsigned char operator[](int i) const {
    u64_t *v=value;
    while (i>=32) {
      i-=32;
      v++;
    }
    return ((*v)>>(i*2)) & 0x3;
  }

  inline operator u64_t() const {
    return *value;
  }

  inline operator unsigned char() const {
    return (*value) & 0x3;
  }

  // hash function for unordered_map
  struct hasher
  {
    inline size_t operator()(binstr_t b) const { return (size_t)*(b.value); }
  };
};

#endif
