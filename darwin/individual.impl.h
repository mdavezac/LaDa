//
//  Version: $Id$
//
#ifndef _INDIVIDUAL_IMPL_H_
#define _INDIVIDUAL_IMPL_H_

#include <opt/debug.h>

namespace Individual
{
  template<class T_INDIVTRAITS>
  void Base<T_INDIVTRAITS> :: clone( const t_This &_indiv )
  {
    quantity   = _indiv.quantity;
    object     = _indiv.object;
    repFitness = _indiv.repFitness;
  }

  template<class T_INDIVTRAITS>
  bool Base<T_INDIVTRAITS> :: operator==( const t_This &_indiv ) const
  {
    if ( Fuzzy::neq(get_concentration(), _indiv.get_concentration() ) )
      return false;
    if ( invalid() or _indiv.invalid() )
      return object == _indiv.object; 
    if ( repFitness != _indiv.repFitness ) return false;
    return object == _indiv.object; 
  }
        
  template<class T_INDIVTRAITS>
  void Base<T_INDIVTRAITS> :: printOn(std::ostream &_os) const
  {
    { // EO stuff
      if (invalid()) {
          _os << "INVALID ";
      }
      else
      {
          _os << repFitness << ' ';
      }
    } // EO stuff
  }

  template<class T_INDIVTRAITS>  template<class SaveOp>
  bool Base<T_INDIVTRAITS> :: Save( TiXmlElement &_node, SaveOp &_saveop ) const
  {
    TiXmlElement *xmlindiv = new TiXmlElement("Individual");
    __DOASSERT( not xmlindiv, "Memory Allocation Error\n" )

    if (     repFitness.Save(*xmlindiv)
         and _saveop(*this, *xmlindiv)  )
    {
      _node.LinkEndChild( xmlindiv );
      return true;
    }

    delete xmlindiv;
    __THROW_ERROR( "Errow while save individual\n" <<  *this << "\n" )
    return false;
  }


  template<class T_INDIVTRAITS>  template<class LoadOp>
  bool Base<T_INDIVTRAITS> ::  Load( const TiXmlElement &_node, LoadOp &_loadop ) 
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if ( name.compare("Individual") )
      parent = _node.FirstChildElement("Individual");
    if ( not parent ) return false;

    if ( not repFitness.Load(*parent) ) return false;

    return _loadop(*this,*parent);
  }

  template<class T_INDIVTRAITS>  template<class SaveOp>
  bool Multi<T_INDIVTRAITS> :: Save( TiXmlElement &_node, SaveOp &_saveop ) const
  {
    TiXmlElement *xmlindiv = new TiXmlElement("Individual");
    __DOASSERT( not xmlindiv, "Memory Allocation Error.\n")

    if (     repFitness.Save(*xmlindiv)
         and _saveop(*this, *xmlindiv)  )
    {
      _node.LinkEndChild( xmlindiv );
      return true;
    }

    delete xmlindiv;
    __THROW_ERROR( "Errow while save individual\n" <<  *this << "\n" )
    return false;
  }

  template<class T_INDIVTRAITS>  template<class LoadOp>
  bool Multi<T_INDIVTRAITS> :: Load( const TiXmlElement &_node, LoadOp &_loadop ) 
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if ( name.compare("Individual") )
      parent = _node.FirstChildElement("Individual");
    if ( not parent ) return false;

    if ( not repFitness.Load(*parent) ) return false;

    return _loadop(*this,*parent);
  }
  
#ifdef _MPI
  template<class T_INDIVTRAITS>  
  bool Base<T_INDIVTRAITS> :: serialize( mpi::BroadCast &_bc )
  {
    return      _bc.serialize( age ) \
            and repFitness.serialize(_bc)\
            and _bc.serialize<t_Object>( object )\
            and t_QuantityTraits::serialize( quantity, _bc );
  }
  template<class T_INDIVTRAITS>  
  bool Base<T_INDIVTRAITS> :: serialize( mpi::BroadCast &_bc ) const
  {
    return      _bc.serialize( age ) \
            and repFitness.serialize(_bc)\
            and _bc.serialize<t_Object>( object )\
            and t_QuantityTraits::serialize( quantity, _bc );
  }
#endif

} // namespace Individual






#endif // _INDIVIDUAL_IMPL_H_
