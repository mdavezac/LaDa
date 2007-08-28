//
//  Version: $Id$
//
#ifndef _ANALYZE_CODE_H_
#define _ANALYZE_CODE_H_

  #ifdef _ANALYZE
    #include <time.h>
    #include <string>
    #include <list>
    #include <fstream>

  
    struct Stamp
    {
      clock_t time;
      std::string name;
 
      Stamp() {};
      Stamp(const Stamp &_stamp ) 
        { time = _stamp.time; name = _stamp.name; };
      Stamp(const std::string  &_name ) 
        { name = _name; time = -clock(); };
      ~Stamp() {};
 
      void operator= (const Stamp &_stamp)
        { name = _stamp.name; time = _stamp.time; } 
      void operator += (const clock_t _time )
        { time += _time; }
      void end()
        { time += clock(); }
      void start () 
        { time = -clock(); }
      void operator() () 
        { time = -clock(); }
      void print_out( ofstream &stream ) const
        { stream << name << " -> " << time << std::endl; }
    };
 
    
    #ifdef CODE_ANALYSIS_IS_MAIN
      std::list<Stamp> stamps;
    #else 
      extern std::list<Stamp> stamps;
    #endif

    #define END_CODE_ANALYSIS(_analysis_file_) \
      std::ofstream analysis_file(_analysis_file_); \
      std::list<Stamp> :: iterator i_stamp = stamps.begin(); \
      std::list<Stamp> :: iterator i_stamp_end = stamps.end(); \
      for( ; i_stamp != i_stamp_end; ++i_stamp ) \
        i_stamp->print_out(analysis_file); \
      analysis_file.flush(); \
      analysis_file.close();

    #define START_ANALYSIS(_codename_) \
      Stamp this_stamp(_codename_);
    #define END_ANALYSIS \
      this_stamp.end(); \
      stamps.push_back( this_stamp );

  #else

    #define END_CODE_ANALYSIS(_analysis_file_)
    #define START_ANALYSIS(_codename_) 
    #define END_ANALYSIS

  #endif
#endif
