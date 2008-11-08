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

namespace LaDa
{
  // Forard declaration.
  //! \cond
  namespace Fitting { class LeaveManyOut; }
  namespace CE
  {
    template< class T_FIT, class T_SOLVER >
      opt::t_ErrorPair leave_many_out( Fitting::LeaveManyOut &_lmo,
                                       const T_FIT &_fit,
                                       const T_SOLVER &_solver );
    namespace Method
    {
      template< class T_COLLAPSE, class T_FIT, class T_MINIMIZER >
        opt::t_ErrorPair leave_many_out( Fitting::LeaveManyOut &_lmo,
                                         T_COLLAPSE &_collapse,
                                         T_FIT &_fit,
                                         const T_MINIMIZER &_min );
    }
  }
  //! \endcond

  namespace Fitting
  {
    //! Creates an object for storing Leave-Many-Out stuff.
    class LeaveManyOut
    {
      template< class T_COLLAPSE, class T_FIT, class T_MINIMIZER >
        friend opt::t_ErrorPair CE::Method::leave_many_out( LeaveManyOut &_lmo,
                                                            T_COLLAPSE &_collapse,
                                                            T_FIT &_fit,
                                                            const T_MINIMIZER &_min );
      template< class T_FIT, class T_SOLVER >
        friend opt::t_ErrorPair CE::leave_many_out( Fitting::LeaveManyOut &_lmo,
                                                    const T_FIT &_fit,
                                                    const T_SOLVER &_solver );
      public:
        //! Type containing an average and a max error.
        typedef std::pair< types::t_real, types::t_real > t_PairErrors;
        //! Type returned by the functor.
        typedef std::pair< opt::ErrorTuple, opt::ErrorTuple > t_Return;

        //! Whether or no to perform.
        bool do_perform;
        //! Level of verbosity.
        types::t_unsigned verbosity;
        
        //! Constructor.
        LeaveManyOut() : do_perform(false), verbosity( 0 ), cmdl_except("except"),
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
        std::vector< types::t_unsigned > except;
        //! The target values to leave out.
        std::vector< std::vector< types::t_unsigned > > sets;
        //! The except command-line string.
        std::string cmdl_except;
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
        ( cmdl_except.c_str(), po::value<std::string>(),
          "A set of structures which should always be in fit (e.g. \" 0 1 2 10 \")" )
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

          if( verbosity >= vlevel1 ) std::cout << "Training:\n";
          intermediate = _interface.check_training( _collapse.function,
                                                    verbosity >= vlevel2 );
          training += opt::ErrorTuple( intermediate.get<1>(), 1e1 );
          if( verbosity >= vlevel1 ) std::cout << intermediate << "\n";

          if( verbosity >= vlevel1 ) std::cout << "Prediction:\n";
          intermediate = _interface.check_predictions( _collapse.function,
                                                       verbosity >= vlevel2 );
          training += opt::ErrorTuple( intermediate.get<1>(), 1e1 );
          if( verbosity >= vlevel1 ) std::cout << intermediate << "\n";
        }

        return t_Return( training, prediction);
      }
      __CATCHCODE(, "Error while performing leave-many-out.\n" )
    }

  }
} // namespace LaDa
#endif
