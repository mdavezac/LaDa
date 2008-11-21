//
//  Version: $Id$
//

#ifndef _LADA_GA_POPULATE_
#define _LADA_GA_POPULATE_

namespace LaDa
{
  namespace GA
  {
    namespace Operators
    {
      //! \brief Grows a population to a set size \a _popsize using operator \a
      //!        _operator to create individuals.
      template< class T_OPERATOR, class T_TABOOS, class T_POPULATION >
        void populate( const T_OPERATOR& _operator,
                       T_TABOOS* const _taboos, T_POPULATION& _pop, 
                       size_t _popsize, size_t _maxiter )
        {
          typedef typename T_POPULATION :: value_type t_Individual;
          types::t_unsigned i = 0, j = 0;
          while ( _pop.size() < _popsize and i < _maxiter)
          {
            t_Individual indiv;
            _operator( indiv );
            if ( _taboos and ( not (*_taboos)( indiv ) ) )
            {
              _pop.push_back(indiv);
                ++j;
            }
            ++i;
          } // while ( i_pop->size() < target )
          __DOASSERT( j < _popsize,
                         "Created " << j << " individuals in " << i << " iterations.\n"
                      << "Are taboos/concentration constraints to restrictive?\n" )
          _pop.resize( _popsize );
        }

      //! Could use a lambda... but can't get it right...
      template< class T_INDIVIDUAL >
        void call_object_mask( T_INDIVIDUAL& _indiv, size_t _first, size_t _end )
          { _indiv.Object().mask( _first, _end ); }

      
      //! \brief Grows a population: (i) by creating an indivdual using \a _initop
      //!        and (ii) adding "masked" individual using \a _maskop.
      template< class T_INITOP, class T_MASKOP, class T_TABOOS, class T_POPULATION >
        void partition_populate( const T_INITOP& _initop, const T_MASKOP& _maskop,
                                 T_TABOOS* const _taboos, T_POPULATION& _pop, 
                                 size_t _popsize, size_t _maxiter )
        {
          typedef typename T_POPULATION :: value_type t_Individual;
          size_t i = 0, popstart( _pop.size() );
          while ( _pop.size() < _popsize and i < _maxiter )
          {
            ++i;
            // first random object
            t_Individual indiv;
            _initop( indiv ); 
            types::t_unsigned start = 0;  
            types::t_unsigned end = indiv.Object().size();  
          
            // indiv: ----
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
            
            _maskop( indiv, start, end );         // indiv: ****
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
            
            _maskop( indiv, start, end/2 );       // indiv: --**
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: **--
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end/4 );       // indiv: -*--
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: *-**
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, end/2, end*3/4 );     // indiv: *--*
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: -**-
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, end/4, end/2 );       // indiv: --*-
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: **-*
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end/2 );       // indiv: ---*
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
            if ( _pop.size() >= _popsize ) break;
       
            _maskop( indiv, start, end );         // indiv: *---
            if ( not ( _taboos and (*_taboos)(indiv) ) ) _pop.push_back(indiv);
          }
          __DOASSERT( _pop.size() < _popsize,
                         "Created " << _pop.size() - popstart << " individuals in "
                      << i << " iterations.\n"
                         "Are taboos/concentration constraints to restrictive?\n" )
          _pop.resize( _popsize );
        }
    } // namespace Operators
  } // namespace GA
} // namespace LaDa

#endif
