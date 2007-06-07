#include <iostream>
#include "arraylist.h"
#include "stringo.h"

class PolytopeFace {
public:
  Array<int> pts;
  Array<Real> up;
  Real c;
};

class LinearInequality {
  public:
   Array<Real> v;
   Real c;
  LinearInequality(int size=0): v(size) {c=0;}
  LinearInequality(const LinearInequality &a): v(a.v) {c=a.c;}
  LinearInequality(const Array<Real> &_v, Real _c): v(_v) {c=_c;}
};

void de_mean(Array<Array<Real> > *px, const Array<Array<Real> > &org_x);

void calc_convex_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x, const Array<Real> &ground);

void update_normals(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &org_x);

int flag_outside_hull(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, int flag_plane_too);

template<class T>
void paste_row_vector(Array<Array<T> > *pout, const Array<Array<T> > &x, const Array<Array<T> > &y);

int is_point_in_hull(const Array<Real> &x, const LinkedList<LinearInequality> &ineq_list);

void clip_hull(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const LinkedList<LinearInequality> &ineq);

void calc_formation(Array<Real> *pfe, Array<Real> *ppure, const Array<Array<Real> > &x, const Array<Real> &e);

void read_inequalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s);

// dim n+1 versions;

template<class T>
void paste_col_vector(Array<Array<T> > *pout, const Array<Array<T> > &x, const Array<T> &y);

void calc_convex_hull_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e);

void update_normals_p1(LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e);

int flag_outside_hull_p1(Array<int> *poutside, LinkedList<PolytopeFace> *hull, const Array<Array<Real> > &x, const Array<Real> &e, int flag_plane_too);
