//
//  Version: $Id$
//
#ifndef __PRINT_OPERATIONS_H_
#define __PRINT_OPERATIONS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iomanip> 

namespace LaDa
{
  namespace Print
  {
    //! \brief Lists standard formatting operations which take no argument.
    //! \details The exact behaviors of the following commands may depend upon
    //!          the class.
    enum t_Operation
    { 
      ENDL,        //!< Generally, will end a line and flush
      FLUSH,       //!< Generally, only flushes
      LEFT,        //!< Generally, left justifies next field.
      RIGHT,       //!< Generally, right justifies next field.
      FIXED,       //!< Generally, sets fixed decimal output for floats and doubles.
      SCIENTIFIC,  //!< Generally, sets scientific decimal output for floats and doubles.
      BOOLALPHA,   //!< Generally, sets bool to ouput as "true" or "false".
      NOBOOLALPHA  //!< Generally, sets bool to ouput as "1" or "0".
    };
    //! Names a standard operation.
    const t_Operation endl        = ENDL;
    //! Names a standard operation.
    const t_Operation flush       = FLUSH;
    //! Names a standard operation.
    const t_Operation left        = LEFT;
    //! Names a standard operation.
    const t_Operation right       = RIGHT;
    //! Names a standard operation.
    const t_Operation fixed       = FIXED;
    //! Names a standard operation.
    const t_Operation scientific  = SCIENTIFIC;
    //! Names a standard operation.
    const t_Operation boolalpha   = BOOLALPHA;
    //! Names a standard operation.
    const t_Operation noboolalpha = NOBOOLALPHA;

    //! \brief Applies standard formatting operations which take no arguments to \a _type.
    //! \details Specialize your template function, as in
    //!          StdOut::operator_<Print::t_Operationt> ( const Print::t_operation& )
    template <class T_TYPE> inline void apply_ops( T_TYPE &_type, const t_Operation &_op)
    {
      switch( _op )
      {
        case ENDL:        _type << std::endl;        break;
        case FLUSH:       _type << std::flush;       break;
        case LEFT:        _type << std::left;        break;
        case RIGHT:       _type << std::right;       break;
        case FIXED:       _type << std::fixed;       break;
        case SCIENTIFIC:  _type << std::scientific;  break;
        case BOOLALPHA:   _type << std::boolalpha;   break;
        case NOBOOLALPHA: _type << std::noboolalpha; break;
      }
    }

    //! \brief Next field will have set width.
    //! \details Specialize your template function, as in
    //!          StdOut::operator_<Print::setw> ( const Print::setw& )
    struct setw
    {
      public:
        types::t_unsigned width; //!< Width of the next field.
      public:
        //! Constructor.
        setw(types::t_unsigned _w) : width(_w) {} 
        //! Stream formatting functor.
        template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
          { _whatever << std::setw(width); }
    };
    //! Next field will have set width.
    typedef setw width;
    //! \brief Next field will have set precision.
    //! \details Specialize your template function, as in
    //!          StdOut::operator_<Print::setprecision> ( const Print::setprecision& )
    class setprecision
    {
      public:
        //! Precision of the next field.
        types::t_unsigned prec;
      public:
        //! Constructor.
        setprecision(types::t_unsigned _n) : prec(_n) {} 
        //! Stream formatting functor.
        template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
          { _whatever << std::setprecision(prec); }
    };
    //! \brief Next field will have set fill character.
    //! \details Specialize your template function, as in
    //!          StdOut::operator_<Print::setfill> ( const Print::setfill& )
    class setfill
    {
      public:
        //! Fill character.
        char c;
      public:
        //! Constructor.
        setfill(char _c) : c(_c) {} 
        //! Stream formatting functor.
        template<class T_TYPE> void operator()( T_TYPE &_whatever ) const
          { _whatever << std::setfill(c); }
    };
  }
} // namespace LaDa
#endif
