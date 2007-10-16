//
//  Version: $Id$
//
#ifndef _STORE_IMPL_H_
#define _STORE_IMPL_H_

#include <print/xmg.h>

namespace Store
{

  template<class T_CONDITION, class T_GATRAITS>
  bool Conditional<T_CONDITION, T_GATRAITS> :: Load (const TiXmlElement &_node)
  {
    if ( not _node.Attribute("print") ) return true;

    std::string str = _node.Attribute("print");
    if( str.find("all") != std::string::npos)
       print_what |= PRINT_RESULTS;
    else if( str.find("results") != std::string::npos ) 
       print_what = PRINT_RESULTS;

    if( str.find("condition") != std::string::npos ) 
       print_what |= PRINT_CONDITION;

    return true;
  }

  template<class T_CONDITION, class T_GATRAITS>
  void Conditional<T_CONDITION, T_GATRAITS> :: operator()( const t_Individual &_indiv )
  {
    // if first evaluation, adds it directly
    if( condition( _indiv ) ) return;

    // checks wether result should be stored
    typename t_Container :: iterator i_found;
#ifdef _MPI
    if ( not new_optima.empty() )
    {
      i_found = std::find( new_optima.begin(), new_optima.end(), _indiv );
      if ( i_found != new_optima.end() ) return;
    }
#endif
    if ( not results.empty() )
    {
      i_found = std::find( results.begin(), results.end(), _indiv );
      if ( i_found != results.end() ) return;
    }
    new_results = true; 
#ifdef _MPI
    if ( not new_optima.empty() ) new_optima.remove_if( condition );
    new_optima.push_back( _indiv );
#else
    if ( not results.empty() ) results.remove_if( condition );
    results.push_back( _indiv );
#endif
    // remove_if should have changed Objective::current_indiv 
    condition( _indiv );
  }

  template<class T_CONDITION, class T_GATRAITS>
  bool Conditional<T_CONDITION, T_GATRAITS> :: Restart( const TiXmlElement &_node )
  {
    const TiXmlElement *xmlresults = &_node;
    std::string name = _node.Value();
    if ( name.compare("Results") ) xmlresults = _node.FirstChildElement("Results");
    if ( not xmlresults ) return false;
    GA::LoadObject<t_GATraits> loadop( evaluator, &t_Evaluator::Load, GA::LOADSAVE_LONG);
    if ( not condition.Restart( *xmlresults, loadop ) )
    {
      Print::xmg << Print::Xmg::comment << "Could not load condition" << Print::endl
                 << Print::Xmg::comment << "Aborting Load of Result"  << Print::endl;
      return false;
    }

    results.clear();
    GA::LoadIndividuals( *xmlresults, loadop, results );
    // This particular line acts as Restart for condition if necessary
    std::remove_if(results.begin(), results.end(), condition); 
    Print::xmg << Print::Xmg::comment << "Reloaded Optimum and "
               << results.size() << " target results" << Print::endl;
    print_results(0, true);
    return true;
  }
  template<class T_CONDITION, class T_GATRAITS>
  bool Conditional<T_CONDITION, T_GATRAITS> :: Save( TiXmlElement &_node ) const
  {
    GA::SaveObject<t_GATraits> saveop( evaluator, &t_Evaluator::Save, GA::LOADSAVE_LONG);
    TiXmlElement *parent = new TiXmlElement("Results");
    if ( not parent ) 
    {
      std::cerr << "Memory error while saving results" << std::endl;
      return false;
    }

    if (      condition.Save(*parent, saveop) 
         and  SaveIndividuals( *parent, saveop, results.begin(), results.end() ) ) 
    {
      _node.LinkEndChild(parent);
      return true;
    }

    std::cerr << "Could not save resuls as requested" << std::endl;
    delete parent;
    return false;
  }

  // Prints either results stored in this object, 
  // or only whatever objective itself does print
  template<class T_CONDITION, class T_GATRAITS>
  void Conditional<T_CONDITION, T_GATRAITS> :: print_results(types::t_unsigned _age,
                                                              bool _is_comment ) const
  {
    t_Base :: print_results( _age, _is_comment );
    if (print_what & PRINT_CONDITION)  
      Print::xmg << ( _is_comment ? Print::make_commented_string( condition.print() ) :
                                    condition.print() )
                 << Print::endl;

    if( not ( print_what & PRINT_RESULTS ) ) return;

    typename t_Container :: const_iterator i_indiv = results.begin();
    typename t_Container :: const_iterator i_end = results.end();
    for (; i_indiv != i_end; ++i_indiv )
      Print::xmg << ( _is_comment ? Print::Xmg::comment: Print::Xmg::clear )
                 << std::setw(12) << std::setprecision(7)
                 << _age << " "
                 << i_indiv->get_concentration() << " "
                 << i_indiv->fitness() << Print::endl;
  }
  template<class T_CONDITION, class T_GATRAITS>
  inline std::string Conditional<T_CONDITION, T_GATRAITS> ::  what_is() const
  { 
    std::ostringstream sstr; sstr << " Conditional on"  << condition.what_is(); 
    return sstr.str();
  }
#ifdef _MPI
  template<class T_CONDITION, class T_GATRAITS>
  bool Conditional<T_CONDITION, T_GATRAITS> :: broadcast( mpi::BroadCast &_bc )
  {
    if ( not condition.broadcast(_bc) ) return false;
    types::t_int n = results.size();
    if ( not _bc.serialize(n) ) return false;
    if ( _bc.get_stage() == mpi::BroadCast::COPYING_FROM_HERE )
      results.resize(n);
    typename t_Container :: iterator i_res = results.begin();
    typename t_Container :: iterator i_res_end = results.end();
    for(; i_res != i_res_end; ++i_res )
      if( not i_res->broadcast( _bc ) ) return false;
    return true;
  }
  template<class T_CONDITION, class T_GATRAITS>
  void Conditional<T_CONDITION, T_GATRAITS> :: synchronize()
  {
    t_Base::synchronize();
    mpi::AllGather allgather( mpi::main );
    typename t_Container :: iterator i_res = new_optima.begin();
    typename t_Container :: iterator i_res_end = new_optima.end();
    for(; i_res != i_res_end; ++i_res )
      i_res->broadcast( allgather );
    allgather.allocate_buffers();
    i_res = new_optima.begin();
    for(; i_res != i_res_end; ++i_res )
      i_res->broadcast( allgather );
    allgather();
    new_optima.clear();
    t_Individual indiv;
    while( indiv.broadcast( allgather ) )
    {
      if( condition( indiv ) ) continue;
      
      // checks wether result should be stored
      typename t_Container :: iterator i_found;
      i_found = std::find( results.begin(), results.end(), indiv );
      if ( i_found != results.end() ) continue;
      results.push_back( indiv );
    }
    results.remove_if( condition );
  }
#endif // _MPI


  namespace Condition
  {
    template< class T_GATRAITS >
    bool BaseOptima<T_GATRAITS> :: Restart( const TiXmlElement &_node, t_LoadOp & _op)
    {
      const TiXmlElement *child = _node.FirstChildElement("optimum");
      if( not child ) return false;
      return optimum.Load( *child, _op );
    }

    template< class T_GATRAITS >
    bool BaseOptima<T_GATRAITS> :: Save( TiXmlElement &_node, t_SaveOp & _op) const
    {
      if ( optimum.invalid() )
      {
        std::cerr << "Optimum is invalid!! when trying to save" << std::endl;
        Print::out << "Optimum is invalid!! when trying to save\n";
        return false;
      }

      TiXmlElement *child = new TiXmlElement("optimum");
      if( not child )
      { 
        std::cerr << "Memory allocation error while save results" << std::endl;
        return false;
      }
      if ( not optimum.Save( *child, _op ) )
      {
        std::cerr << "Could not save optimum" << std::endl;
        return false;
      }

      _node.LinkEndChild(child);
      return true;
    }

    template< class T_GATRAITS >
    std::string BaseOptima<T_GATRAITS> :: print() const
    {
      std::ostringstream sstr;
      sstr << optimum << std::setw(12) << std::setprecision(7) 
           << optimum.fitness();
      return sstr.str(); 
    }

    template< class T_GATRAITS >
    FromObjective<T_GATRAITS> :: FromObjective   ( const TiXmlElement &_node ) 
                                               : t_Base(_node), objective(NULL),
                                                 val(0), end_val(0), delta(0), owns_objective(false)
    {
      double d=0;
      if ( _node.Attribute("delta") )
        _node.Attribute("delta", &d);
      delta = (t_ScalarQuantity) std::abs(d);

      objective = t_Objective :: new_from_xml( _node );
      if ( not objective )
        throw std::runtime_error("Could not Load objective from input in conditional store\n"); 
      owns_objective = true;
    }
    template< class T_GATRAITS >
    FromObjective<T_GATRAITS> :: FromObjective   ( typename t_Objective::Vector* _type,
                                                   const TiXmlElement &_node )
                                               : t_Base(_node), objective(_type),
                                                 val(0), end_val(0), delta(0),
                                                 owns_objective(false), print_condition(false)
    {
      double d=0;
      if ( _node.Attribute("delta") )
        _node.Attribute("delta", &d);
      delta = (t_ScalarQuantity) std::abs(d);

      if( _node.Attribute("print") )
      {
        std::string name = _node.Attribute("print");
        if (    name.compare("all") == 0
             or name.compare("condition") == 0 ) print_condition=true;
      }
 
      if ( not objective )
        throw std::runtime_error("Could not Load objective from input in conditional store\n"); 
    }

    template< class T_GATRAITS >
    bool FromObjective<T_GATRAITS> :: operator()( const t_Individual &_indiv )
    {
      objective->init( _indiv );
      if ( optimum.invalid() )
      {
        optimum = _indiv;
        val = (*objective)( optimum.const_quantities() );
        end_val = val + delta;
        return false;
      }
      t_ScalarQuantity indiv_val = (*objective)(_indiv.const_quantities());
      if ( t_QuantityTraits::less(indiv_val, val) )
      {
        optimum = _indiv;
        val = indiv_val;
        end_val = val + delta;
        return false;
      }

      return t_QuantityTraits::greater(indiv_val, end_val); 
    }
    
    template< class T_GATRAITS >
    inline std::string FromObjective<T_GATRAITS> :: what_is() const
    { 
      std::ostringstream sstr;
      sstr << objective->what_is() << " Objective, with delta = " << delta;
      return sstr.str();
    } 

    template< class T_GATRAITS >
    bool FromObjective<T_GATRAITS> :: Save( TiXmlElement &_node, t_SaveOp & _op) const
    {
      if (     t_Base::Save( _node, _op)  
           and objective->Save( _node, _op) )  return true;
      
      std::cerr << "Could not save objective" << std::endl;
      return false;
    }
  } // namespace Condition

} // namespace Store

#endif // _STORE_IMPL_H_
