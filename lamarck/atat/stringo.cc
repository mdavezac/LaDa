#include "stringo.h"


namespace atat
{ 

char * AutoString::empty="";
/*
istream& operator >> (istream &file, AutoString &str) {
	static char buf[1000];
	file.get(buf,1000);
	return file;
}
*/

Real to_real(const AutoString &s) {
  istringstream str((const char *)s);
  Real r;
  str >> r;
  return r;


} // namespace atat
}
