#ifndef __STRINGO_H__
#define __STRINGO_H__

#include <string.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "misc.h"

namespace atat
{ 


class AutoString {
protected:
	char *buf;
	static char *empty;
public:
	AutoString(const char *str) {
	  buf=NULL;
	  set(str);
	}
	AutoString(const AutoString &a) {
	  buf=NULL;
	  set(a);
	}
	AutoString(int len) {
	  buf=NULL;
	  set(len);
	}
	void set(const char *str) {
	  if (buf) delete[] buf;
	  if (str) {
	    buf=new char[strlen(str)+1];
	    strcpy(buf,str);
	  }
	  else
	    buf=NULL;
	}
	void set(int len) {
	  if (buf) delete[] buf;
	  buf=new char[len+1];
	  MEMCLR(buf,len+1);
	}
	AutoString(void) {
	  buf=NULL;
	}
	int len(void) const {
	  if (buf) {
	    return strlen(buf);
	  }
	  else {
	    return 0;
	  }
	}
	void operator=(const AutoString &a) {
	  set(a);
	}
#ifndef STRING_FIX
	operator char * () {
	  if (buf)
	    return buf;
	  else
	    return empty;
	}
#endif
	operator const char * () const {
	  if (buf)
	    return buf;
	  else
	    return empty;
	}
	char& operator [](int i) {
#ifdef DEBUG
	  if ( i<0 || i>=len() ) {
	    cerr << "AutoString out of range: " << i << "/" << len() << endl;
	    return *empty;
	  }
	  else
#endif
	    return buf[i];
	}
	const char& operator [](int i) const {
#ifdef DEBUG
	  if ( i<0 || i>=len() ) {
	    cerr << "AutoString out of range: " << i << "/" << len() << endl;
	    return *empty;
	  }
	  else
#endif
	    return buf[i];
	}
	void operator += (const AutoString &s) {
	  char *newbuf=new char[len()+s.len()+1];
	  strcpy(newbuf,*this);
	  strcpy(newbuf+len(),s);
	  delete[] buf;
	  buf=newbuf;
	}
	void operator += (char c) {
	  char *newbuf=new char[len()+2];
	  strcpy(newbuf,*this);
	  newbuf[len()]=c;
	  newbuf[len()+1]=0;
	  delete[] buf;
	  buf=newbuf;
	}
	int operator == (const AutoString &s) const {
	  return (strcmp(*this,s)==0);
	}
	~AutoString() {
	  if (buf) delete[] buf;
	}
	friend istream& operator >> (istream &file, AutoString &str);
};

Real to_real(const AutoString &s);


} // namespace atat

#endif
