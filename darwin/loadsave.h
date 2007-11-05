//
//  Version: $Id$
//
#ifndef _DARWIN_LOADSAVE_H_
#define _DARWIN_LOADSAVE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

namespace GA
{

  //! \brief Wrapper around GA::Evaluator derived member %function for saving
  //!        individuals from XML.
  //! \details Saving an individual, and specifically an object is obviously
  //!          application-specifics. Furthermore, we want to be able to save
  //!          an individual both in a short condensed form when saving is for
  //!          internal %GA reasons, and in a long declarative form which can
  //!          be understood by a user. As such, Individual::Base takes a
  //!          functor as argument when saving itself. This is is the functor.
  //!          It wraps itself around an  GA::Evaluator member %function which
  //!          does the actual work.
  template<class T_GATRAITS>
  class SaveObject
  {
    public:
      //! All %types relevant to %GA.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! Type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! Type of the application-specific evaluator class
      typedef typename t_GATraits :: t_Evaluator  t_Evaluator;
      //! Type of the member routine which will do the saving in the evaluator class 
      typedef bool(t_Evaluator::*t_saveop)( const t_Individual &_indiv,
                                            TiXmlElement &_node, bool _type ) const;

    public:
      //! Instance to an application-specific evaluator class
      t_Evaluator &evaluator;      
      //! Pointer to the member %function which will do the saving.
      t_saveop op;      
      ///! Wether to save in long or condensed form
      bool type;

      //! \brief Constructor and Initializer.
      //! \param _e Reference to the evaluator
      //! \param _op Pointer to the member %function of \a _e which will do the saving.
      //! \param _t Whether to save in long or short form.
      SaveObject   ( t_Evaluator &_e, t_saveop _op, bool _t)
                 : evaluator(_e), op(_op), type(_t) {}
      //! Functor routine
      bool operator()(const t_Individual &_indiv, TiXmlElement &_node ) const
        { return (evaluator.*op)(_indiv, _node, type); }
  };

  //! \brief Wrapper around GA::Evaluator derived member %function for loading
  //!        individuals from XML. 
  //! \details Loading an individual, and specifically an object is obviously
  //!          application-specifics. Furthermore, we want to be able to load
  //!          an individual both in a short condensed form when loading is for
  //!          internal %GA reasons, and in a long declarative form which can
  //!          be understood by a user. As such, Individual::Base takes a
  //!          functor as argument when loading itself. This is is the functor.
  //!          It wraps itself around an  GA::Evaluator member %function which
  //!          does the actual work.
  template<class T_GATRAITS>
  class LoadObject
  {
    public:
      //! All %types relevant to %GA.
      typedef T_GATRAITS t_GATraits;
    protected:
      //! Type of the individual
      typedef typename t_GATraits :: t_Individual t_Individual;
      //! Type of the application-specific evaluator class
      typedef typename t_GATraits :: t_Evaluator  t_Evaluator;
      //! Type of the member routine which will do the loading in the evaluator class 
      typedef bool( t_Evaluator::*t_loadop )( t_Individual &_indiv,
                                              const TiXmlElement &_node, bool _type );

    public:
      //! Instance to an application-specific evaluator class
      t_Evaluator &evaluator;      
      //! Pointer to the member %function which will do the loading.
      t_loadop op;      
      ///! Wether to load in long or condensed form
      bool type;

      //! \brief Constructor and Initializer.
      //! \param _e Reference to the evaluator
      //! \param _op Pointer to the member %function of \a _e which will do the loading.
      //! \param _t Whether to load in long or short form.
      LoadObject   ( t_Evaluator &_e, t_loadop _op, bool _t)
                 : evaluator(_e), op(_op), type(_t) {}
      //! Functor routine
      bool operator()(t_Individual &_indiv, const TiXmlElement &_node )
        { return (evaluator.*op)(_indiv, _node, type); }
  };

  //! Alias for saving/loading in short condensed form.
  const bool LOADSAVE_SHORT = true;
  //! Alias for saving/loading in long declarative form.
  const bool LOADSAVE_LONG = false;

  //! \brief Helper %function for saving a whole range of individuals.
  //! \param _node Node below which to save the range of individuals.
  //! \param _saveop Functor for saving an individual.
  //! \param _first First individual in the range
  //! \param _end one past the last individual in the range
  //! \pre The range [ \a _first, \a _end ) should be valid.
  template<class T_IT, class SaveOp>
  bool SaveIndividuals( TiXmlElement &_node, SaveOp &_saveop, T_IT _first, T_IT _end) 
  {
    bool result = true;
    for(; _first != _end; ++_first )
      if( not _first->Save( _node, _saveop ) ) result = false;
    return result;
  }


  //! \brief Helper %function for loading a whole container of individuals.
  //! \details In the default version (\a _n=0), this %function keeps on
  //!          loading individuals until it gets an error or cannot find
  //!          another individual tag within \a _node. If on the other hand \a
  //!          _n is not null, then the routines exits after loading the first
  //!          \a _n individuals or when it finds no more individuals if there
  //!          are less than \a _n.
  //! \param _node Node from which to load the individuals.
  //! \param _loadop Functor for loading an individual.
  //! \param _container Functor for loading an individual.
  //! \param _n the maximum number of individuals to load. Zero means load all
  //1           you can find. Zero is the default.
  template<class T_CONTAINER, class LoadOp>
  bool LoadIndividuals( const TiXmlElement &_node, LoadOp &_loadop, 
                        T_CONTAINER& _container, types::t_unsigned _n = 0 )
  {
    bool did_load = false;
    const TiXmlElement *indiv_xml = &_node;
    std::string name = _node.Value();
    if( name.compare("Individual") )
      indiv_xml = _node.FirstChildElement("Individual");
    if( not indiv_xml )
      return false;
      
    for(; indiv_xml; indiv_xml = indiv_xml->NextSiblingElement("Individual") )
    {
      if ( _n and _container.size() > _n )  break;
      typename T_CONTAINER :: value_type indiv;
      if ( indiv.Load( *indiv_xml, _loadop ) )
      {
        _container.push_back( indiv );
        did_load = true; 
      }
    }
    return did_load;
  }

} // namespace GA
#endif // _DARWIN_LOADSAVE_H_
