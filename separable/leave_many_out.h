//
//  Version: $Id$
//
#ifndef _LEAVE_MANY_OUT_H_
#define _LEAVE_MANY_OUT_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <string>
#include <sstream>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include <opt/types.h>
#include <opt/debug.h>
#include <opt/errors.h>

namespace Fitting
{
  class LeaveManyOut
  {
    public:
      //! Type containing an average and a max error.
      typedef std::pair< types::t_real, types::t_real > t_PairErrors;
      //! Type returned by the functor.
      typedef std::pair< opt::ErrorTuple, opt::ErrorTuple > t_Return;

      //! Whether or no to perform.
      bool do_perform;
      //! Level of verbosity.
      types::t_unsigned verbose;
      
      //! Constructor.
      LeaveManyOut() : do_perform(false), verbose( false ), 
                       cmdl_set( "sets" ), cmdl_file("setsname"),
                       filename( "lmo_sets" ), nb_sets(0), nb_perset(0) {}
      //! Destructor..
      ~LeaveManyOut() {}


      //! Load from command-line.
      void add_cmdl( boost::program_options::options_description &_options );
      //! Extract command-line options.
      void extract_cmdl( boost::program_options::variables_map &_map );

      //! Functor.
      template< class T_INTERFACE, class T_ALLSQ, class T_COLLAPSE>
      t_Return operator()( T_INTERFACE &_interface,
                           T_ALLSQ &_allsq, 
                           T_COLLAPSE &_collapse );

    protected:

      //! Creates a set of sets.
      void create_sets( types::t_unsigned _tsize );

      //! Reads set of sets from input.
      void read_sets();

      //! The target values to leave out.
      std::vector< std::vector< types::t_unsigned > > sets;
      //! The set command-line string.
      std::string cmdl_set;
      //! The input filename command-line string.
      std::string cmdl_file;
      //! The file where to save/load a set.
      std::string filename;
      //! Number of sets.
      types::t_unsigned nb_sets;
      //! Number per set.
      types::t_int nb_perset;
      //! \cond
      const bool static vlevel1 = 1; 
      const bool static vlevel2 = 3; 
      //! \endcond
  };

  inline void LeaveManyOut :: add_cmdl( boost::program_options
                                             ::options_description &_options )
  {
    __DEBUGTRYBEGIN
    namespace po = boost::program_options;
    _options.add_options()
      ( cmdl_set.c_str(), po::value<std::string>(),
        "Defines the size and number "
        "of leave-many-out sets \" NUMBER-SETS  NUMBER-PER-SET \"." )
      ( cmdl_file.c_str(), po::value<std::string>(),
        "Name of a file where to save/load the leave-many-out sets." );
    __DEBUGTRYEND(, "Error while creating LeaveManyOut command-line options.\n" );
  }

  template< class T_INTERFACE, class T_ALLSQ, class T_COLLAPSE>
  LeaveManyOut :: t_Return LeaveManyOut :: operator()( T_INTERFACE &_interface,
                                                       T_ALLSQ &_allsq, 
                                                       T_COLLAPSE &_collapse )
  {
    opt::ErrorTuple training;
    opt::ErrorTuple prediction;
    if( not do_perform ) return t_Return( training, prediction );
    try
    {
      if( not sets.size() ) create_sets( _interface.nb_structs() );
      _interface.exclude.resize(1, 0);

      bool first_iter = true;
      typedef std::vector< std::vector< types::t_unsigned > >
                                    :: const_iterator const_iterator;
      const_iterator i_set = sets.begin();
      const_iterator i_set_end = sets.end();
      opt::ErrorTuple error;
      for(; i_set != i_set_end; ++i_set, first_iter=false )
      {
        opt::ErrorTuple intermediate;
        _interface.exclude.resize( i_set->size() );
        std::copy( i_set->begin(), i_set->end(), _interface.exclude.begin() );

        _interface.fit( _allsq, _collapse );

        if( verbose >= vlevel1 ) std::cout << "Training:\n";
        intermediate = _interface.check_training( _collapse.function,
                                                  verbose >= vlevel2 );
        training += opt::ErrorTuple( intermediate.get<1>(), 1e1 );
        if( verbose >= vlevel1 ) std::cout << intermediate << "\n";

        if( verbose >= vlevel1 ) std::cout << "Prediction:\n";
        intermediate = _interface.check_predictions( _collapse.function,
                                                     verbose >= vlevel2 );
        training += opt::ErrorTuple( intermediate.get<1>(), 1e1 );
        if( verbose >= vlevel1 ) std::cout << intermediate << "\n";
      }

      return t_Return( training, prediction);
    }
    __CATCHCODE(, "Error while performing leave-many-out.\n" )
  }

}

#endif
