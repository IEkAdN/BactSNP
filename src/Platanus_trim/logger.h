// -*- mode: c++ -*-
#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <iostream>
#include <stdarg.h>
#include "binstr.h"
using namespace std;

/**
 * stream wrapper for information and debug messages.
 *
 * usage:
 *   string str = "string";
 *   Logger logger;
 *   logger << "hoge=" << str << endl;
 *   logger.debug() << "hoge=" << str << endl;
 *   logger.error().printf("hoge=0x%x\n",255);
 */
class Logger {
protected:
  enum LEVEL {DEBUG=0, INFO, NOTICE, ERROR};
  ostream&  ostr;
  LEVEL     level, current;
  unsigned  format;

public:
  Logger(LEVEL l=NOTICE,ostream& o=cerr):ostr(o),level(l),current(l),format(0) {}

  template<typename T>
  Logger& operator<<(T t) { if(level<=current)ostr<<t; return *this; }

  template<typename T>
  Logger& operator<<(const T *t) { if(level<=current)ostr<<t; return *this; }

  Logger& operator<<(u64_t& t) { return (*this)<<(u64_t)t; }
  Logger& operator<<(u64_t t) {
    int i;
    switch (this->format) {
    case 0: // hex
      this->printf("0x%016lx", t);
      break;
    case 1: // bin
      for (i=62; i>=0; i-=2) {
        ostr << (int)((t>>i)&0x3);
      }
      break;
    case 2: // char
      for (i=62; i>=0; i-=2) {
        ostr << "ACGT"[(t>>i)&0x3];
      }
    }
    return *this;
  }
  Logger& operator<<(binstr_t& t) { return (*this)<<(binstr_t)t; }
  Logger& operator<<(binstr_t t)  { 
    if (level<=current) {
      int i, nb=(t.len+31)/32;
      u64_t *v=t.value+nb-1;
      char buf[] = "0x%0..lx";
      switch (this->format) {
      case 0:  // hex
        snprintf(buf, sizeof(buf), "0x%%0%dlx", 1+(t.len%32)/2);
        this->printf(static_cast<const char*>(buf), *v);
        while (v--!=t.value) {
          this->printf("|%016lx", *v);
        }
        break;
      case 1:  // bin
        for (i=(t.len%32-1)*2; i>=0; i-=2) {
          ostr << (int)(((*v)>>i)&0x3);
        }
        while (v--!=t.value) {
          ostr << "|";
          for (i=62; i>=0; i-=2) {
            ostr << (int)(((*v)>>i)&0x3);
          }
        }
        break;
      case 2:  // char (N->0)
        for (i=(t.len%32-1)*2; i>=0; i-=2) {
          ostr << "ACGT"[((*v)>>i)&0x3];
        }
        while (v--!=t.value) {
          ostr << "|";
          for (i=62; i>=0; i-=2) {
            ostr << "ACGT"[((*v)>>i)&0x3];
          }
        }
      }
    }
    return *this;
  }
  Logger& operator<<(string& s)   { if(level<=current)ostr<<s.c_str(); return *this; }
  Logger& operator<<(ostream& manip(ostream&)) { if(level<=current)manip(ostr); return *this; }

  Logger& printf(const char* fmt, ...) {
    va_list args;
    va_start(args, fmt);
    if (level<=current) {
      char buf[1024];
      vsnprintf(buf, sizeof(buf), fmt, args);
      ostr << buf;
    }
    va_end(args);
    return *this;
  }


  // for change level
  Logger& debug() { current=DEBUG; return *this; }
  Logger& info() { current=INFO; return *this; }
  Logger& notice() { current=NOTICE; return *this; }
  Logger& error() { current=ERROR; return *this; }
  Logger& reset() { current=level; return *this; }

  // for change format
  Logger& hex() { format=0; return *this; }
  Logger& bin() { format=1; return *this; }
  Logger& chr() { format=2; return *this; }
};

#endif
