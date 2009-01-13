//
//  Version: $Id: checkpoints.h 849 2008-11-12 04:45:34Z davezac $
//
#ifndef _LADA_GA_CHECKPOINTS_PRINT_POPULATION_H_
#define _LADA_GA_CHECKPOINTS_PRINT_POPULATION_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <print/xmg.h>
#include <print/stdout.h>
#include <opt/modifiers.h>
#include <opt/tinyxml.h>

#include "../gencount.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! \brief Prints any population container. 
      //! \details Can be a container of rvalues, pointers, or reference to individuals.
      struct PrintPopContainer
      {
        //! Title of the container.
        const std::string title;
        //! Constructor.
        PrintPopContainer( const std::string& _title ) : title( _title ) {}
        //! Copy Constructor.
        PrintPopContainer( const PrintPopContainer& _c ) : title( _c.title ) {}

        //! Functor. 
        template< class T_CONTAINER >
          void operator()( bool _lastcall, const T_CONTAINER& _container )
          {
            if( _lastcall ) return;
            Print::out << title << "\n";
            typename T_CONTAINER :: const_iterator i_indiv = _container.begin();
            typename T_CONTAINER :: const_iterator i_indiv_end = _container.end();
            for(; i_indiv != i_indiv_end; ++i_indiv )
              Print::out << Modifier::const_innermost( *i_indiv ) << "  Fitness: " 
                         << Print::fixed << Print::setw(12) << Print::setprecision(5) << "  "
                         << Modifier::const_innermost(*i_indiv).fitness() << "\n";
            Print::out << "\n";
          }
      };
      //! \brief Prints new individuals.
      //! \details Can be a container of rvalues, pointers, or reference to individuals.
      struct PrintNewIndividuals
      {
        //! Reference to the generation counter.
        const GenCount age;
        //! Constructor.
        PrintNewIndividuals( const GenCount _age ) : age( _age ) {}
        //! Copy Constructor.
        PrintNewIndividuals( const PrintNewIndividuals& _c ) : age( _c.age ) {}

        //! Functor. 
        template< class T_CONTAINER >
          void operator()( bool _lastcall, const T_CONTAINER& _container )
          {
            if( _lastcall ) return;
            if ( age() == 0 ) return;
            Print::out << "New Individuals:\n";
            typename T_CONTAINER :: const_iterator i_indiv = _container.begin();
            typename T_CONTAINER :: const_iterator i_indiv_end = _container.end();
            for(; i_indiv != i_indiv_end; ++i_indiv )
              if( age() == Modifier::const_innermost( *i_indiv ).get_age() )
                Print::out << Modifier::const_innermost( *i_indiv ) << "  Fitness: " 
                           << Print::fixed << Print::setw(12) << Print::setprecision(5) << "  "
                           << Modifier::const_innermost(*i_indiv).fitness() << "\n";
            Print::out << "\n";
          }
      };

      namespace Factory
      {
        //! Factory function for printing offspring at end of generation.
        template< class T_CHECKPOINT >
          void print_newindividuals( T_CHECKPOINT& _checkpoint, 
                                     const GenCount _age )
          {
            _checkpoint.connect_statistic( PrintNewIndividuals( _age ) ); 
            Print::out << "Will print new individuals at each generation.\n";
            Print::xmg << Print::Xmg::comment << "Will print new individuals at each generation." 
                       << Print::endl;
          }
        //! Factory function for printing population at end of generation.
        template< class T_CHECKPOINT >
          void print_population( T_CHECKPOINT& _checkpoint )
          {
            _checkpoint.connect_statistic( PrintPopContainer( "Current Population:" ) ); 
            Print::out << "Will print current population at each generation.\n";
            Print::xmg << Print::Xmg::comment << "Will print current population at each generation." 
                       << Print::endl;
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::print_offspring()
        template< class T_CHECKPOINT >
          void (*print_newindividuals( const T_CHECKPOINT& ))( T_CHECKPOINT&, const GenCount ) 
            { return &Factory::print_newindividuals< T_CHECKPOINT >; }
        //! Returns address of void LaDa::GA::CheckPoint::Factory::print_offspring()
        template< class T_CHECKPOINT >
          void (*print_population( const T_CHECKPOINT& ))( T_CHECKPOINT& ) 
            { return &Factory::print_population< T_CHECKPOINT >; }
      }
  
    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif
