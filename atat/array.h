#ifndef __ARRAY_H__
#define __ARRAY_H__

#include <iostream>
#include "vectmac.h"
#include "misc.h"

#include <opt/types.h>

namespace atat
{ 

template<class T>
class Array {
 protected:
  types::t_int size;
  T *buf;
 protected:
  void init(types::t_int new_size) {
    size=new_size;
    if (new_size!=0)
      buf=new T[size];
    else
      buf=NULL;
  }
  T& access(types::t_int i) const {
#ifdef DEBUG
    if (i<0 || i>=size) {
      cerr << "Array out of range: " << i << "/" << size << endl;
      COREDUMP
      return buf[0];
    }
    else
#endif
      return buf[i];
  }
 public:
  Array(types::t_int new_size=0) {
    init(new_size);
  }
  Array(const Array<T> &a) {
    init(a.size);
    for (types::t_int i=0; i<a.size; i++) {
      buf[i]=a.buf[i];
    }
  }
  ~Array(void) {
    delete[] buf;
  }
  types::t_int operator()(void) const {
    return size;
  }
  types::t_int getSize(void) const {
    return size;
  }
  types::t_int get_size(void) const {
    return size;
  }
  T& operator()(types::t_int i) {
    return access(i);
  }
  const T& operator()(types::t_int i) const {
    return access(i);
  }
  T& operator[](types::t_int i) {
    return access(i);
  }
  const T& operator[](types::t_int i) const {
    return access(i);
  }
  void operator=(const Array<T> &a) {
    resize(a.size);
    for (types::t_int i=0; i<a.size; i++) {
      buf[i]=a.buf[i];
    }
  }
  void copy(const Array<T> &a) {
    resize(a.size);
    for (types::t_int i=0; i<a.size; i++) {
      buf[i]=a.buf[i];
    }
  }
  void resize(types::t_int new_size) {
    if (new_size!=size) {
      delete[] buf;
      init(new_size);
    }
  }
  T* get_buf(void) {
    return buf;
  }
  /*
  operator T* () {
    return buf;
  }
  operator const T* () const {
    return buf;
  }
  */
};

template<class T>
void zero_array(Array<T> *a) {
  for (types::t_int i=0; i<a->getSize(); i++) {
    (*a)(i)=(T)0;
  }
}

template<class T>
void one_array(Array<T> *a) {
  for (types::t_int i=0; i<a->getSize(); i++) {
    (*a)(i)=(T)1;
  }
}

template<class T>
types::t_int operator==(const Array<T> &a, const Array<T> &b) {
  if (a.get_size()!=b.get_size()) return 0;
  for (types::t_int i=0; i<a.get_size(); i++) {
    if (a(i)!=b(i)) return 0;
  }
  return 1;
}

template<class T>
types::t_int operator!=(const Array<T> &a, const Array<T> &b) {
  if (a.get_size()!=b.get_size()) return 1;
  for (types::t_int i=0; i<a.get_size(); i++) {
    if (a(i)!=b(i)) return 1;
  }
  return 0;
}

class Array2dIterator;

template<class T>
class Array2d {
 protected:
  iVector2d size;
  T *buf;
 protected:
  void init(const iVector2d& new_size) {
    size=new_size;
    types::t_int surface=size(0)*size(1);
    if (surface>0) {
      buf=new T[surface];
    }
    else
      buf=NULL;
  }
  T& access(types::t_int x, types::t_int y) const {
#ifdef DEBUG
    if (x<0 || x>=size(0) || y<0 || y>=size(1)) {
      cerr << "Array2d out of range: " << x << " " << y << " / " << size(0) << " " << size(1) << endl;
      COREDUMP
      return buf[0];
    }
    else
#endif
      return buf[y*size(0)+x];
  }
 public:
  Array2d(void) {
    init(iVector2d(0,0));
  }
  Array2d(const iVector2d& new_size) {
    init(new_size);
  }
  Array2d(types::t_int new_size_x, types::t_int new_size_y) {
    init(iVector2d(new_size_x,new_size_y));
  }
  Array2d(const Array2d<T> &a) {
    init(a.get_size());
    for (types::t_int i=0; i<size(0)*size(1); i++) {
      buf[i]=a.buf[i];
    }
  }
  void operator=(const Array2d<T> &a) {
    init(a.get_size());
    for (types::t_int i=0; i<size(0)*size(1); i++) {
      buf[i]=a.buf[i];
    }
  }
  ~Array2d(void) {
    delete[] buf;
  }
  iVector2d operator()(void) const {
    return size;
  }
  iVector2d getSize(void) const {
    return size;
  }
  iVector2d get_size(void) const {
    return size;
  }
  T& operator()(iVector2d i) {
    return access(i(0),i(1));
  }
  T& operator()(types::t_int x, types::t_int y) {
    return access(x,y);
  }
  const T& operator()(iVector2d i) const {
    return access(i(0),i(1));
  }
  const T& operator()(types::t_int x, types::t_int y) const {
    return access(x,y);
  }
  /*  void copy(const Array2d<T> &a) {
    resize(a.getSize());
    Array2dIterator i(a.getSize());
    while (i) {
      THIS(i)=a(i);
      i++;
    }
  } */
  void resize(iVector2d new_size) {
    if (new_size!=size) {
      delete[] buf;
      init(new_size);
    }
  }
  void resize(types::t_int new_size_x, types::t_int new_size_y) {
    resize(iVector2d(new_size_x,new_size_y));
  }
  operator T* () {
    return buf;
  }
  operator const T* () const {
    return buf;
  }
};

class Array2dIterator {
  iVector2d i;
  iVector2d size;
public:
  Array2dIterator(iVector2d _size) {
    init(_size);
  }
  void init(iVector2d _size) {
    i=iVector2d(0,0);
    size=_size;
  }
  void init(void) {
    i=iVector2d(0,0);
  }
  types::t_int operator++(int) {
    i(0)++;
    if (i(0)==size(0)) {
      i(0)=0;
      i(1)++;
      if (i(1)==size(1)) {
	//i(1)=0;
	return 0;
      }
    }
    return 1;
  }
  types::t_int operator+=(iVector2d step) {
    i(0)+=step(0);
    if (i(0)>=size(0)) {
      i(0)=0;
      i(1)+=step(1);
      if (i(1)>=size(1)) {
	//i(1)=0;
	return 0;
      }
    }
    return 1;
  }
  operator iVector2d(void) const {
    return i;
  }
  operator types::t_int(void) const {
    if (i(1)>=size(1) || i(0)>=size(0))
      return 0;
    else
      return 1;
  }
  types::t_int eol(void) {return (i(0)==(size(0)-1));}
};

template<class T>
void zero_array(Array2d<T> *a) {
  Array2dIterator i(a->getSize());
  while (i) {
    (*a)(i)=(T)0;
    i++;
  }
}

template<class T>
ostream& operator << (ostream &s, const Array<T> &a) {
  s << a.get_size() << endl;
  for (types::t_int i=0; i<a.get_size(); i++) {
    s << a(i) << endl;
  }
  return s;
}

template<class T>
istream& operator >> (istream &s, Array<T> &a) {
  types::t_int size;
  s >> size;
  a.resize(size);
  for (types::t_int i=0; i<a.get_size(); i++) {
    s >> a(i);
  }
  return s;
}

template<class T>
ostream& operator << (ostream &s, const Array2d<T> &a) {
  iVector2d size=a.get_size();
  s << size << endl;
  for (types::t_int i=0; i<size(0); i++) {
    for (types::t_int j=0; j<size(1); j++) {
      s << a(i,j) << " ";
    }
    s << endl;
  }
  return s;
}

template<class T>
istream& operator >> (istream &s, Array2d<T> &a) {
  iVector2d size;
  s >> size;
  a.resize(size);
  for (types::t_int i=0; i<size(0); i++) {
    for (types::t_int j=0; j<size(1); j++) {
      s >> a(i,j);
    }
  }
  return s;
}

template<class T>
types::t_int find_value(const Array<T> &a, const T &x) {
  types::t_int minj=-1;
  types::t_real mind=MAXFLOAT;
  for (types::t_int j=0; j<a.get_size(); j++) {
    types::t_real d=fabs((types::t_real)(a(j)-x));
    if (d<mind) {
      mind=d;
      minj=j;
    }
  }
  return minj;
}

template<class T>
inline T max(const Array<T> &a) {
  T m=a(0);
  for (types::t_int i=1; i<a.get_size(); i++) {
    if (a(i)>m) m=a(i);
  }
  return m;
}

template<class T>
inline T min(const Array<T> &a) {
  T m=a(0);
  for (types::t_int i=1; i<a.get_size(); i++) {
    if (a(i)<m) m=a(i);
  }
  return m;
}

template<class T>
inline types::t_int index_max(const Array<T> &a) {
  T m=a(0);
  types::t_int idx=0;
  for (types::t_int i=1; i<a.get_size(); i++) {
    if (a(i)>m) {m=a(i); idx=i;}
  }
  return idx;
}

template<class T>
inline types::t_int index_min(const Array<T> &a) {
  T m=a(0);
  types::t_int idx=0;
  for (types::t_int i=1; i<a.get_size(); i++) {
    if (a(i)<m) {m=a(i); idx=i;}
  }
  return idx;
}

template<class T>
types::t_int is_in_array(const Array<T> &a, const T &x) {
  for (types::t_int i=0; i<a.get_size(); i++) {
    if (x==a(i)) return 1;
  }
  return 0;
}

template<class T>
types::t_int index_in_array(const Array<T> &a, const T &x) {
  for (types::t_int i=0; i<a.get_size(); i++) {
    if ((T)x==(T)(a(i))) return i;
  }
  return -1;
}

template<class T>
void sort_array(Array<T> *a) {
  for (types::t_int i=0; i<a->get_size()-1; i++) {
    for (types::t_int j=i; j<a->get_size()-1; j++) {
      if ((*a)(j)>(*a)(j+1)) {swap(&(*a)(j),&(*a)(j+1));}
    }
  }
}

template<class T>
void resize(Array<Array<T> > *a, types::t_int n1, types::t_int n2){
  a->resize(n1);
  for (types::t_int i1=0; i1<n1; i1++) {
    (*a)(i1).resize(n2);
  }
}

template<class T>
void ArrayArray_to_Array2d(Array2d<T> *pa2d, const Array<Array<T> > &a, types::t_int permute) {
  if (permute) {
    pa2d->resize(a(0).get_size(),a.get_size());
    for (types::t_int i=0; i<pa2d->get_size()(0); i++) {
      for (types::t_int j=0; j<pa2d->get_size()(1); j++) {
	(*pa2d)(j,i)=a(i)(j);
      }
    }
  }
  else {
    pa2d->resize(a.get_size(),a(0).get_size());
    for (types::t_int i=0; i<pa2d->get_size()(0); i++) {
      for (types::t_int j=0; j<pa2d->get_size()(1); j++) {
	(*pa2d)(i,j)=a(i)(j);
      }
    }
  }
}

typedef Array<types::t_real> ArrayReal;
typedef Array<types::t_int> Arrayint;
typedef Array<Arrayint> ArrayArrayint;
typedef Array<ArrayReal> ArrayArrayReal;

template<class T>
class ArrayFunctionArray {
 public:
  virtual void eval(Array<T> *py, const Array<T> &x) {}
  virtual types::t_int read(istream &file) {return 0;}
};


} // namespace atat

#endif
