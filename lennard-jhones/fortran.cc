//
//  Version: $Id$
//
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <fstream>

#include <boost/spirit/include/classic_primitives.hpp>
#include <boost/spirit/include/classic_numerics.hpp>
#include <boost/spirit/include/classic_parser.hpp>
#include <boost/spirit/include/classic_actions.hpp>
#include <boost/spirit/include/classic_kleene_star.hpp>
#include <boost/spirit/include/classic_kleene_star.hpp>
#include <boost/spirit/include/classic_optional.hpp>
#include <boost/spirit/include/classic_assign_actor.hpp>
#include <boost/spirit/include/classic_sequence.hpp>
#include <boost/spirit/include/classic_operators.hpp>
#include <boost/bind.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "fortran.h"
#include "c2f_handles.h"


namespace LaDa
{
  namespace CLJ
  {
    void Relaxer :: read_info( const std::string& _filename )
    {
      __TRYBEGIN
      namespace fs = boost::filesystem;
      namespace bsc = boost::spirit::classic;


      std::string name = "gsl_bfgs2";
      types::t_real tolerance(1e-6);
      types::t_real linetolerance(1e-2);
      types::t_real linestep(1e-1);
      size_t itermax(50);
      bool verbose(false);
      std::string strategy("fast");

      std::ifstream ldas( _filename.c_str(), std::ifstream::in );
      std::string line;
      while( std::getline( ldas, line ) )
      {
        if( line.find( "#" ) != std::string::npos )
          line = line.substr( 0, line.find( "#" ) );
        bsc::parse
        (
          line.c_str(),
          bsc::str_p("type:") >> !(*bsc::space_p)
                              >> (*(bsc::alnum_p | bsc::str_p('_')))[ bsc::assign_a( name ) ]
        );
        bsc::parse
        (
          line.c_str(),
          bsc::str_p("tolerance:") - bsc::str_p("line tolerance")
            >> !(*bsc::space_p) >> bsc::real_p[ bsc::assign_a( tolerance ) ]
        );
        bsc::parse
        (
          line.c_str(),
          bsc::str_p("line tolerance:")
            >> !(*bsc::space_p) >> bsc::real_p[ bsc::assign_a( linetolerance ) ]
        );
        bsc::parse
        (
          line.c_str(),
          bsc::str_p("line step:")
            >> !(*bsc::space_p) >> bsc::real_p[ bsc::assign_a( linestep ) ]
        );
        bsc::parse
        (
          line.c_str(),
          bsc::str_p("itermax:")
            >> !(*bsc::space_p) >> bsc::real_p[ bsc::assign_a( itermax ) ]
        );
        if
        (
          bsc::parse
          (
            line.c_str(),
            bsc::str_p("verbose")
          ).hit
        ) verbose = true;
        bsc::parse
        (
          line.c_str(),
          "strategy:" >> !(*bsc::space_p) >> (*bsc::alpha_p)[ bsc::assign_a( strategy ) ]
        );
      }

      TiXmlElement fakexml( "Minimizer" );
      fakexml.SetAttribute( "type", name );
      fakexml.SetDoubleAttribute( "tolerance", tolerance );
      fakexml.SetAttribute( "itermax", itermax );
      fakexml.SetDoubleAttribute( "linetolerance", linetolerance );
      fakexml.SetDoubleAttribute( "linestep", linestep );
      fakexml.SetAttribute( "strategy", strategy );
      fakexml.SetAttribute( "verbose", verbose ? "true": "false" );
      __DOASSERT( not minimizer_.Load( fakexml ), 
                  "Could not create minimizer.\n" << fakexml << "\n" )

      __TRYEND(,"Error while parsing " << _filename << "\n" )
    }


    Relaxer::t_Type Relaxer :: operator()( const size_t _Natoms,
                                           const int *const _occupation, 
                                           t_Type *const _cell, 
                                           t_Type *const _positions, 
                                           t_Type *const _stress, 
                                           t_Type *const _forces ) 
    {
       const int natoms( _Natoms );
       t_Wrapper :: t_Arg args( _Natoms );
       Binder_ bound( functional_, natoms, _occupation );
       wrapper_.init( bound );
       wrapper_.create_arguments( args, _Natoms, _cell, _positions );
       minimizer_( wrapper_, args );
       return wrapper_.results( args, _cell, _positions, _stress, _forces );
    }

  } // namespace CLJ.
} // namespace LaDa

extern "C" void FC_FUNC_(create_relaxer, CREATE_RELAXER)
           ( 
             int* _handle,
             const int* _fsize,
             const char* _filename,
             LaDa :: CLJ :: Relaxer :: t_Functional _function
           )
{
  __TRYBEGIN
    __DOASSERT( *_fsize <= 0, "Filename-length negative or zero." )
    std::string in;
    std::copy( _filename, _filename + *_fsize, std::back_inserter( in ) );

    LaDa::CLJ::Relaxer relaxer;
    relaxer.read_info( in );
    relaxer.init( _function );
    LaDa::C2FHandles< LaDa::CLJ::Relaxer > handles;
    *_handle = handles.push_back( relaxer );
  __TRYEND(, "Could not create relaxer." )
}

extern "C" void FC_FUNC_(release_relaxer, RELEASE_RELAXER)( int* const _handle )
{
  __TRYBEGIN
    __DOASSERT( *_handle < 0, "Handle cannot be negative.\n" )
    LaDa::C2FHandles< LaDa::CLJ::Relaxer > handles;
    handles.erase( *_handle );
  __TRYEND(, "Could not delete relaxer." )
}

extern "C" void FC_FUNC_(call_relaxer, CALL_RELAXER)
           ( 
             int* const _handle,             // Relaxer handle
             const int* const _natoms,       // Natoms
             double* const    _cell,         // Cell
             const int* const _occupations,  // occupations
             double* const    _positions,    // positions
             double* const    _forces,       // forces
             double* const    _stress,       // stress
             double*          _energy        // energy
           )
{
  __DOASSERT( *_handle < 0, "Handle cannot be negative.\n" )
  LaDa::C2FHandles< LaDa::CLJ::Relaxer > handles;
  LaDa::CLJ::Relaxer& relaxer( handles[ *_handle ] ); 
  *_energy = relaxer( *_natoms, _occupations, _cell, _positions, _stress, _forces );
}
extern "C" double FC_FUNC_(boost_erfc, BOOST_ERFC)( const double *const _in )
  { return boost::math::erfc( *_in ); }
