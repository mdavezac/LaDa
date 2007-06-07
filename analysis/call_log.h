#ifndef _CALL_LOG_H_
#define _CALL_LOG_H_

  #ifdef CALL_LOG
    #include <string>
    #include <ext/hash_map>
    #include <fstream>
    #include <functional>
    #include <opt/types.h>
    using __gnu_cxx::hash_map;
    using __gnu_cxx::hash;

    struct compare_strings
    {
      bool operator()(const string &s1, const string &s2) const
      {
        return s1.compare(s2) == 0;
      }
      size_t operator()(const string &s1) const
      {
        hash<const char*> H;
        return H( s1.c_str() );
      }
    };
  
    struct Call_Log
    {
      hash_map<string, types::t_int, compare_strings> stamps;
 
      Call_Log() { stamps.clear(); };
      ~Call_Log() {};
 
      void operator() ( const char *stamp_name )
      { ( stamps[stamp_name] )++; }
      void operator() ( std::string &stamp_name )
      { ( stamps[stamp_name] )++; }
      void print_log( const char *_analysis_file_ )
      {
        std::ofstream analysis_file(_analysis_file_); 
        hash_map<string, types::t_int, compare_strings> :: iterator i_stamp = stamps.begin(); 
        hash_map<string, types::t_int, compare_strings> :: iterator i_stamp_end = stamps.end(); 
        for( ; i_stamp != i_stamp_end; ++i_stamp ) 
        {
          analysis_file << "procedure : " << i_stamp->first 
                        << " calls: "  << i_stamp->second << std::endl;
        }
        analysis_file.flush(); 
        analysis_file.close();
      }
      types::t_int value(std::string &stamp_name)
      {
        hash_map<string, types::t_int, compare_strings> :: iterator i_stamp = stamps.find( stamp_name ); 
        return i_stamp->second;
      }

    };
 
    
    #ifdef CODE_LOG_IS_MAIN
      Call_Log call_log;
    #else
      extern Call_Log call_log;
    #endif

    #define PRINT_LOG(_analysis_file_) \
      call_log.print_log(_analysis_file_); 
    #define LOG_PROCEDURE(_codename_) \
      call_log(_codename_);

  #else
    #define PRINT_LOG(_analysis_file_) 
    #define LOG_PROCEDURE(_codename_) 

  #endif
#endif
