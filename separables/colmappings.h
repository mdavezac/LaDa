#ifndef _CE_COLMAPPINGS_H_
#define _CE_COLMAPPINGS_H_

#include "LaDaConfig.h"

#include<boost/lambda/bind.hpp>
#include<boost/numeric/ublas/vector.hpp>
#include<vector>

#include <misc/types.h>
#include <opt/debug.h>

namespace LaDa
{
  namespace CE
  {
    namespace Mapping
    {
      class Basic
      {
        public:
          //! Type of vectors.
          typedef std::vector<types::t_real> t_Vector;

          //! Constructor.
          Basic() {}
          //! Copy Constructor.
          Basic  (const Basic& _c)
                : N(_c.N), weights_( _c.weights_ ),
                  targets_( _c.targets_ ) {}
          //! Destructor
          ~Basic() {}
          //! Returns weights.
          t_Vector& weights() { return weights_; }
          //! Returns weights.
          const t_Vector& weights() const { return weights_; }
          //! Returns one weight.
          t_Vector::value_type weight( size_t i ) const
          { 
            LADA_NASSERT( N != weights_.size(), "Incoherent logic.\n" )
            LADA_NASSERT( i >= weights_.size(), 
                         "Index out-of-range: " << i
                      << " >= " << weights_.size() << "\n" )
            return weights_[i]; 
          }
          //! Returns weights.
          t_Vector& targets() { return targets_; }
          //! Returns weights.
          const t_Vector &targets() const { return targets_; }
          //! Returns one weight.
          t_Vector::value_type target( size_t i ) const 
          {
            LADA_NASSERT( N != targets_.size(), "Incoherent logic.\n" )
            LADA_NASSERT( i >= targets_.size(), 
                         "Index out-of-range: " << i
                      << " >= " << targets_.size() << "\n" )
            return targets_[i]; 
          }
          //! Returns the number of structures.
          size_t size() const { return N; }
          //! Allows to skip out on a structure for leave-one or many-out.
          bool do_skip( size_t _i ) const { return false; }
          //! Initializes the mapping.
          template< class T_STRUCTURES >
            void init( const T_STRUCTURES& _strs );

        protected:
          //! Number of structures.
          size_t N;
          //! Weights of structures.
          t_Vector weights_;
          //! Target values of structures.
          t_Vector targets_;
      };

      class SymEquiv : public Basic
      {
        public:
          //! Type of vectors.
          typedef Basic::t_Vector t_Vector;

          //! Constructor.
          SymEquiv() {}
          //! Copy Constructor.
          SymEquiv  (const Basic& _c)
                  : Basic( _c ) {} 
          //! Copy Constructor.
          SymEquiv  (const SymEquiv& _c)
                  : Basic( _c ),
                    equiweights( _c.equiweights ), 
                    nb_( _c.nb_ ) {} 
          //! Destructor
          ~SymEquiv() {}
          //! Returns the range of equivalent configurations for structure \a _i.
          boost::numeric::ublas::range range( size_t _i ) const 
          {
            LADA_NASSERT( _i >= nb_.size(), "Index out-of-range.\n" )
            LADA_NASSERT( _i + 1 >= nb_.size(), "Index out-of-range.\n" )
            return boost::numeric::ublas::range( nb_[_i], nb_[_i+1] ); 
          }
          //! Returns the weight for an equivalent structure.
          t_Vector::value_type eweight( size_t _i, size_t _c ) const
          {
            LADA_NASSERT( _i + 1 >= nb_.size(), "Index out-of-range.\n" )
            LADA_NASSERT( _c >= nb_[_i+1] - nb_[_i], "Index out-of-range.\n" )
            LADA_NASSERT( _i + _c >= equiweights.size(), "Index out-of-range.\n" )
            return equiweights[ nb_[_i] + _c ]; 
          }
          //! Initializes the mapping.
          template< class T_STRUCTURES, class T_CONFIGURATIONS >
            void init( const T_STRUCTURES& _strs, 
                       const T_CONFIGURATIONS& _confs );


        protected:
          //! Weights of structures.
          t_Vector equiweights;
          //! Structure ranges.
          std::vector< size_t > nb_;
      };

    template< class T_BASE, class T_CONTAINER = std::vector< types::t_unsigned>  >
      class ExcludeMany : public T_BASE
      {
        public:
          //! Type of the base class.
          typedef T_BASE t_Base;
          //! Type of the set of excluded.
          typedef T_CONTAINER t_Container;
          //! Index of structure to exclude.
          T_CONTAINER excluded;
          //! Constructor.
          ExcludeMany() : T_BASE() {}
          //! Copy Constructor.
          ExcludeMany( const T_BASE & _c ) : T_BASE( _c ) {}
          //! Copy Constructor.
          ExcludeMany   ( const ExcludeMany & _c )
                      : T_BASE( _c ), excluded( _c.excluded ) {}
          //! Destructor.
          ~ExcludeMany() {}
          //! Returns true if \a _i in ExcludeMany::excluded.
          bool do_skip( size_t _i ) const
          { 
            if( excluded.empty() ) return false;
            return std::find( excluded.begin(), excluded.end(), _i ) != excluded.end();
          }
      };

    template< class T_BASE, class T_CONTAINER >
      class ExcludeMany<T_BASE, T_CONTAINER*> : public T_BASE
      {
        public:
          //! Type of the base class.
          typedef T_BASE t_Base;
          //! Type of the set of excluded.
          typedef T_CONTAINER* t_Container;
          //! Index of structure to exclude.
          T_CONTAINER *excluded;
          //! Constructor.
          ExcludeMany() : T_BASE(), excluded(NULL) {}
          //! Copy Constructor.
          ExcludeMany( const T_BASE & _c ) : T_BASE( _c ), excluded(NULL) {}
          //! Copy Constructor.
          ExcludeMany   ( const ExcludeMany & _c )
                      : T_BASE( _c ), excluded( _c.excluded ) {}
          //! Destructor.
          ~ExcludeMany() {}
          //! Returns true if \a _i in ExcludeMany::excluded.
          bool do_skip( size_t _i ) const
          { 
            if( excluded->empty() ) return false;
            return std::find( excluded->begin(), excluded->end(), _i ) != excluded->end();
          }
      };

    template< class T_BASE >
      class ExcludeOne : public T_BASE
      {
        public:
          //! Type of the base class.
          typedef T_BASE t_Base;
          //! Index of structure to exclude.
          size_t n;
          //! Wether to exclude at all.
          bool do_exclude;
          //! Constructor.
          ExcludeOne() : T_BASE(), n(0), do_exclude( false ) {}
          //! Copy Constructor.
          ExcludeOne( const T_BASE & _c ) : T_BASE( _c ), n(0), do_exclude( false ) {}
          //! Copy Constructor.
          ExcludeOne   ( const ExcludeOne & _c )
                      : T_BASE( _c ), n( _c.n ) {}
          //! Destructor.
          ~ExcludeOne() {}
          //! Returns true if \a _i == ExcludeOne::n and ExcludeOne::do_exclude is true.
          bool do_skip( size_t _i ) const { return do_exclude and _i == n; }
      };


    template< class T_STRUCTURES >
      void Basic :: init( const T_STRUCTURES& _strs )
      {
        namespace bl = boost::lambda;
        // Copy structural weights first.
        weights_.clear();
        std::transform
        (
          _strs.begin(), _strs.end(), std::back_inserter(weights_),
          bl::bind( &T_STRUCTURES :: value_type :: weight, bl::_1 )
        );

        // Copy structural energies second.
        targets_.clear();
        std::transform
        (
          _strs.begin(), _strs.end(), std::back_inserter(targets_),
          bl::bind( &T_STRUCTURES :: value_type :: energy, bl::_1 )
        );
        N = _strs.size();
      }

    template< class T_STRUCTURES, class T_CONFIGURATIONS >
      void SymEquiv::init( const T_STRUCTURES& _strs, 
                           const T_CONFIGURATIONS& _confs )
      {
        namespace bl = boost::lambda;
        
        Basic::init( _strs );

        // Construct internal weights (between equivalent confs) 
        // and initializes the nb_ structure.
        nb_.resize(1, 0);
        equiweights.clear();
        size_t sum(0);
        typename T_CONFIGURATIONS :: const_iterator i_confs = _confs.begin();
        typename T_CONFIGURATIONS :: const_iterator i_confs_end = _confs.end();
        for(; i_confs != i_confs_end; ++i_confs )
        {
          sum += i_confs->size();
          nb_.push_back( sum );
          typename T_CONFIGURATIONS :: value_type 
                                    :: const_iterator i_conf = i_confs->begin();
          typename T_CONFIGURATIONS :: value_type 
                                    :: const_iterator i_conf_end = i_confs->end();
          for(; i_conf != i_conf_end; ++i_conf )
            equiweights.push_back( i_conf->second );
        }
      }

    } // end of Mapping namespace.

  } // end of CE namespace
} // namespace LaDa
#endif
