#ifndef _DARWIN_LOADSAVE_H_
#define _DARWIN_LOADSAVE_H_

#include <tinyxml/tinyxml.h>

namespace darwin
{

  template<class T_OBJECT, class T_EVALUATOR>
  class SaveObject
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_EVALUATOR t_Evaluator;
    private:
      typedef bool( t_Evaluator::*t_saveop )( const t_Object &_object, TiXmlElement &_node, bool _type ) const;

    public:
      t_Evaluator &evaluator;      
      t_saveop op;      
      bool type;

      SaveObject   ( t_Evaluator &_e, t_saveop _op, bool _t)
                 : evaluator(_e), op(_op), type(_t) {}
      bool operator()(const t_Object &_object, TiXmlElement &_node ) const
        { return (evaluator.*op)(_object, _node, type); }
  };
  template<class T_OBJECT, class T_EVALUATOR>
  class LoadObject
  {
    public:
      typedef T_OBJECT t_Object;
      typedef T_EVALUATOR t_Evaluator;
    private:
      typedef bool( t_Evaluator::*t_loadop )( t_Object &_object, const TiXmlElement &_node, bool _type );

    public:
      t_Evaluator &evaluator;      
      t_loadop op;      
      bool type;

      LoadObject   ( t_Evaluator &_e, t_loadop _op, bool _t)
                 : evaluator(_e), op(_op), type(_t) {}
      bool operator()(t_Object &_object, const TiXmlElement &_node )
        { return (evaluator.*op)(_object, _node, type); }
  };

  const bool LOADSAVE_SHORT = true;
  const bool LOADSAVE_LONG = false;

  template<class T_IT, class SaveOp>
  void SaveIndividuals( TiXmlElement &_node, SaveOp &_saveop, T_IT _first, T_IT _end) 
  {
    for(; _first != _end; ++_first )
      _first->Save( _node, _saveop );
  }
  template<class T_CONTAINER, class LoadOp>
  bool LoadIndividuals( const TiXmlElement &_node, LoadOp &_loadop, 
                        T_CONTAINER& _container, types::t_unsigned _size = 0 )
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
      if ( _size and _container.size() > _size )
        break;
      typename T_CONTAINER :: value_type indiv;
      if ( indiv.Load( *indiv_xml, _loadop ) )
      {
        _container.push_back( indiv );
        did_load = true; 
      }
    }
    return did_load;
  }

} // namespace darwin
#endif // _DARWIN_LOADSAVE_H_
