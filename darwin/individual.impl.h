//
//  Version: $Id$
//
#ifndef _INDIVIDUAL_IMPL_H_
#define _INDIVIDUAL_IMPL_H_

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
    if ( std::abs(get_concentration() - _indiv.get_concentration()) > types::tolerance )
      return false;
    if ( invalid() or _indiv.invalid() )
      return object == _indiv.object; 
    if ( repFitness != _indiv.repFitness ) return false;
    return object == _indiv.object; 
  }
        
  template<class T_INDIVTRAITS>
  const typename Base<T_INDIVTRAITS> :: t_Fitness& Base<T_INDIVTRAITS> :: fitness() const
  {
     if ( (bool) invalid() )
       throw std::runtime_error("invalid fitness\n");
     return repFitness;
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
    if ( not xmlindiv )
    {
      std::cerr << "Memory Allocation Error while saving individual" << std::endl;
      return false;
    }

    if (     repFitness.Save(*xmlindiv)
         and _saveop(*this, *xmlindiv)  )
    {
      _node.LinkEndChild( xmlindiv );
      return true;
    }

    delete xmlindiv;
    std::cerr << "Errow while save individual" << std::endl <<  *this << std::endl;
    return false;
  }


  template<class T_INDIVTRAITS>  template<class LoadOp>
  bool Base<T_INDIVTRAITS> ::  Load( const TiXmlElement &_node, LoadOp &_loadop ) 
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if ( name.compare("Individual") )
      parent = _node.FirstChildElement("Individual");
    if ( not parent )
      return false;

    if ( not repFitness.Load(*parent) )
      return false;

    return _loadop(*this,*parent);
  }

#ifdef _MPI
  template<class T_INDIVTRAITS> 
  bool Base<T_INDIVTRAITS> :: broadcast( mpi::BroadCast &_bc )
  {
    return      _bc.serialize( age ) 
            and repFitness.broadcast(_bc)
            and _bc.serialize<t_Object>( object )
            and t_QuantityTraits::broadcast( quantity, _bc );
  }
#endif

  template<class T_INDIVTRAITS>  template<class SaveOp>
  bool Multi<T_INDIVTRAITS> :: Save( TiXmlElement &_node, SaveOp &_saveop ) const
  {
    TiXmlElement *xmlindiv = new TiXmlElement("Individual");
    if ( not xmlindiv )
    {
      std::cerr << "Memory Allocation Error while saving individual" << std::endl;
      return false;
    }

    if (     repFitness.Save(*xmlindiv)
         and _saveop(*this, *xmlindiv)  )
    {
      _node.LinkEndChild( xmlindiv );
      return true;
    }

    delete xmlindiv;
    std::cerr << "Errow while save individual" << std::endl <<  *this << std::endl;
    return false;
  }

  template<class T_INDIVTRAITS>  template<class LoadOp>
  bool Multi<T_INDIVTRAITS> :: Load( const TiXmlElement &_node, LoadOp &_loadop ) 
  {
    const TiXmlElement *parent = &_node;
    std::string name = parent->Value();
    if ( name.compare("Individual") )
      parent = _node.FirstChildElement("Individual");
    if ( not parent )
      return false;

    if ( not repFitness.Load(*parent) )
      return false;

    return _loadop(*this,*parent);
  }
  
} // namespace Individual

#endif // _INDIVIDUAL_IMPL_H_
