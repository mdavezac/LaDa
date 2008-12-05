//
//  Version: $Id$
//
#ifndef _LADA_GA_CHECKPOINTS_SAVE_EVERY_H_
#define _LADA_GA_CHECKPOINTS_SAVE_EVERY_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <boost/signals.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include <tinyxml/tinyxml.h>

#include <print/stdout.h>
#include <print/xmg.h>
#include <mpi/mpi_object.h>
#include <opt/debug.h>
#include "../taboos/history.h"
#include "../loadsave.h"

namespace LaDa
{
  namespace GA
  {
    namespace CheckPoint
    {
      //! Holds all saving functors.
      class SaveGAState
      {
        public:
          //! Constructor.
          SaveGAState() : signals_(new t_SaveSignals) {}
          //! Copy Constructor.
          SaveGAState( const SaveGAState & _c ) : signals_(_c.signals_) {}
          //! Updater functor.
          void operator()( bool );
          //! Connects a save signal.
          template< class T_FUNCTOR > void connect( T_FUNCTOR &_functor )
            { signals_->connect( _functor ); }
           
        protected:
          //! Type of saving sockets.
          typedef boost::signal< void( TiXmlElement& ) > t_SaveSignals;
          //! Save sockets.
          boost::shared_ptr<t_SaveSignals> signals_;
          //! Filename where to save.
          boost::filesystem::path filename_;
      };

      //! Saves a functor.
      template< class T_FUNCTOR > 
        void save_from_functor( TiXmlElement& _node, 
                                const std::string _childname,
                                const T_FUNCTOR _functor )
        {
          __TRYBEGIN
            TiXmlElement *child = new TiXmlElement(_childname);
            _functor( *child );
            _node.LinkEndChild( child );
            Print::out << ("Saved " + _childname + "\n");
          __TRYEND(, ("Could not save " + _childname + "\n") )
        }

      //! Saves a population.
      template< class T_FUNCTOR, class T_CONTAINER >
        void save_population( TiXmlElement& _node,
                              const std::string _childname,
                              const T_FUNCTOR _saveop, 
                              const T_CONTAINER& _container )
        {
          __TRYBEGIN
            save_from_functor
            (
              _node, "Population", 
              boost::bind
              ( 
                &GA::SaveIndividuals< typename T_CONTAINER::const_iterator, T_FUNCTOR >,
                _1, _saveop, _container.begin(), _container.end() 
              )
            );
          __TRYEND(, ("Could not save container " + _childname + "\n") )
        }

      //! Saves all islands.
      template< class T_FUNCTOR, class T_ISLANDS >
        void save_islands( TiXmlElement& _node, 
                           const T_FUNCTOR _saveop,
                           const T_ISLANDS& _islands )
        {
          __TRYBEGIN
            TiXmlElement *child = new TiXmlElement("Population");
            foreach( const typename T_ISLANDS :: value_type pop, _islands )
              save_population( "Island", _saveop, pop );
            _node.LinkEndChild( child );
            Print::out << ("Saved Islands.\n");
          __TRYEND(, ("Could not save Islands.\n") )
        }

      //! Save a store functor.
      template< class T_DARWIN >
        void save_store( const T_DARWIN& _darwin )  { _darwin.store->Save(); }

      namespace Factory
      {

        //! Factory function for stopping on finding file.
        template< class T_CHECKPOINT, class T_DARWIN >
          void save_every( T_CHECKPOINT& _checkpoint,
                           const TiXmlElement& _node, 
                           const T_DARWIN &_darwin  )
          {
            const boost::filesystem::path path( _darwin.save_filename );
            if( not _node.Attribute( "what" ) ) return;
            const std::string what( _node.Attribute( "what" ) );
            bool saveall( what.compare( "all" ) or what.compare( "ALL" ) );
            typedef typename T_DARWIN :: t_Evaluator t_Evaluator;
            typedef typename t_Evaluator :: t_GATraits t_GATraits;
            typedef typename t_GATraits :: t_Individual t_Individual;
            typedef SaveObject< t_GATraits > t_SaveOp;
            t_SaveOp saveop( _darwin.evaluator, &t_Evaluator::Save, LOADSAVE_SHORT);
            SaveGAState signals;
            if( saveall or what.find( "history" ) != std::string::npos )
            {
              signals.connect
              (
                 boost::bind
                 (
                   &GA::CheckPoint::save_population< t_SaveOp, Taboo::History<t_Individual> >,
                   _1, "History", saveop, boost::cref(_darwin.history)
                 )
              );
              Print::out << "Will Save History.\n";
            }
            if( saveall or what.find( "results" ) != std::string::npos )
            {
              signals.connect
              (
                 boost::bind
                 (
                   &GA::CheckPoint::save_store< T_DARWIN >,
                   _1, _darwin
                 )
              );
              Print::out << "Will Save Results.\n";
            }
            if( saveall or what.find( "pop" ) != std::string::npos )
            {
              signals.connect
              (
                 boost::bind
                 (
                   &GA::CheckPoint::save_islands< t_SaveOp, typename T_DARWIN :: t_Islands >,
                   _1, saveop, boost::cref(_darwin.islands)
                 )
              );
              Print::out << "Will Save Results.\n";
            }
            _checkpoint.connect_updater( signals );
          }
      }
      namespace AddressOf
      {
        //! Returns address of void LaDa::GA::CheckPoint::Factory::max_generations()
        template< class T_CHECKPOINT, class T_DARWIN >
          void (*save_every( const T_CHECKPOINT&, const T_DARWIN& ))
              ( T_CHECKPOINT&, const TiXmlElement&, const T_DARWIN& ) 
            { return &Factory::save_every< T_CHECKPOINT, T_DARWIN>; }
      }
  

    } // namespace CheckPoint
  } // namespace GA
} // namespace LaDa

#endif