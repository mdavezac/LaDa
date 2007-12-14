//
//  Version: $Id$
//
#ifndef _PRINT_MANIP_H_
#define _PRINT_MANIP_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <cctype>

//! Contains print out related stuff.
namespace Print
{

  //! Removes everything from the beginning to the last "/" character in the string.
  std::string StripDir( std::string _string );
  //! \brief Removes everything in \a _dir found after the first occurence of
  //!        pattern \a _str..
  //! \details Fisrt reduces \a _dir to the string before the first occurence
  //!          of \a _str. Then removes everything after the first occurrence
  //!          of either "/", "\t", or "\n".
  std::string StripDir( const std::string &_dir, const std::string &_str );
  //! \brief Erase everything before the first occurrence and after the first
  //!        occurrence of any " \t\n" character.
  std::string StripEdges( std::string _string );
  //! \brief Replaces the substring between the first "~" and the next "/" with the
  //!        home directory of the currrent user.
  std::string reformat_home( std::string _str );
  //! Transforms all characters to lowercase.
  std::string lowercase(std::string _string);
  //! Transforms all characters to uppercase.
  std::string uppercase(std::string _string);


  //! Returns true if c is uppercase
  inline bool is_uppercase( char c ) { return c != std::tolower( c ); }
  //! Returns true if c is lowercase
  inline bool is_lowercase( char c ) { return c != std::toupper( c ); }
}
#endif 
