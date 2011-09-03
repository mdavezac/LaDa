#ifndef LADA_COMPARE_SITES_H
#define LADA_COMPARE_SITES_H

#include "LaDaConfig.h"

#include <algorithm>

#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/mpl/or.hpp>
#include <boost/foreach.hpp>

#include "atom.h"
#include "is_container.h"
#include "utilities.h"


namespace LaDa
{
  namespace crystal 
  {
    //! Compares site occupations.
    template<class T_TYPE, class ENABLE=void> struct CompareOccupations;
 
    //! Comparison of scalar types.
    template<class T_TYPE> 
      struct CompareOccupations
        < 
          T_TYPE,
          typename boost::enable_if<details::is_scalar<T_TYPE> >::type
        >
      {
        public:
          //! Constructor.
          CompareOccupations(T_TYPE const &_ref) : reference_(_ref) {}
          //! Copy constructor.
          CompareOccupations(CompareOccupations const &_c) : reference_(_c.reference_) {}
          //! Performs comparison.
          bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
        protected:
          //! Object to which to compare.
          T_TYPE reference_;
          //! Compares string or integer occupation.
          template<class T>
            typename boost::enable_if
              < 
                typename boost::mpl::or_
                  < typename boost::is_convertible<T, std::string> :: type,
                    typename boost::is_integral<T> > :: type,
                bool
              >::type cmp_(T const &_type) const { return reference_ == _type; }
          //! Compares floating point occupations.
          template<class T>
            typename boost::enable_if<boost::is_floating_point<T>, bool>::type
              cmp_(T const &_type) const { return math::eq(_type, reference_); }
      };

    //! Comparison of vector types.
    template<class T_TYPE> 
      struct CompareOccupations
        < 
          T_TYPE,
          typename boost::enable_if<details::is_container<T_TYPE> >::type
        >
      {
        public:
          //! Constructor.
          CompareOccupations(T_TYPE const &_ref) : reference_(_ref) {}
          //! Copy constructor.
          CompareOccupations(CompareOccupations const &_c) : reference_(_c.reference_) {}
          //! Performs comparison.
          bool operator()( T_TYPE const& _type ) const { return cmp_<T_TYPE>(_type); }
        protected:
          //! Object to which to compare.
          T_TYPE reference_;
          //! Compares string or integer occupations.
          template<class T>
            typename boost::enable_if
              < 
                typename boost::mpl::or_
                  < boost::is_convertible<typename T::value_type, std::string>,
                    boost::is_integral<typename T::value_type> > :: type,
                bool
              >::type cmp_( T const &_type ) const
              { 
                foreach(typename T::value_type const &_t, _type)
                  if( std::find(reference_.begin(), reference_.end(), _t) == reference_.end() ) return false;
                foreach(typename T_TYPE::value_type const &_t, reference_)
                  if( std::find(_type.begin(), _type.end(), _t) == _type.end() ) return false;
                return true;
              }
          //! Compares floating point occupations.
          template<class T>
            typename boost::enable_if< boost::is_floating_point<typename T::value_type>, bool>::type
              cmp_( T const _type ) const
              {
                foreach(typename T::value_type const &_t, _type)
                {
                  typename T::const_iterator const i_found
                    = std::find_if( reference_.begin(), reference_.end(),
                                    CompareOccupations<typename T::value_type>(_t) );
                  if( i_found == reference_.end() ) return false;
                }
                foreach(typename T_TYPE::value_type const &_t, reference_)
                {
                  typename T::const_iterator const i_found
                    = std::find_if( _type.begin(), _type.end(),
                                    CompareOccupations<typename T::value_type>(_t));
                  if( i_found == _type.end() ) return false;
                }
                return true;
              }
      };

    //! Comparison of set types.
    template<class T_TYPE> 
      struct CompareOccupations
        < 
          T_TYPE,
          typename boost::enable_if<details::is_set<T_TYPE> >::type
        >
      {
        public:
          //! Constructor.
          CompareOccupations(T_TYPE const &_ref) : reference_(_ref) {}
          //! Copy constructor.
          CompareOccupations(CompareOccupations const &_c) : reference_(_c.reference_) {}
          //! Performs comparison.
          bool operator()( T_TYPE const& _type ) const { return _type == reference_; }
        protected:
          //! Object to which to compare.
          T_TYPE reference_;
      };

    //! Comparison of positions.
    struct ComparePositions
    {
      public:
        //! Constructor.
        template<class T_DERIVED>
        ComparePositions   ( Eigen::DenseBase<T_DERIVED> const &_ref, types::t_real _tol=types::tolerance)
                         : reference_(_ref), tolerance_(_tol) {}
        //! Copy constructor.
        ComparePositions(ComparePositions const &_c) : reference_(_c.reference_), tolerance_(_c.tolerance_) {}
        //! Performs comparison.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const
            { return math::eq(_pos, reference_, tolerance_); }
      protected:
        //! Object to which to compare.
        math::rVector3d reference_;
        //! Tolerance for comparison.
        types::t_real tolerance_;
    };
    
    //! Comparison positions with tolerance a function of distance.
    struct RelComparePositions
    {
      public:
        //! Constructor.
        template<class T_DERIVED>
        RelComparePositions   ( Eigen::DenseBase<T_DERIVED> const &_ref,
                                types::t_real _linear=types::tolerance,
                                types::t_real _constant=types::tolerance )
                            : reference_(_ref), linear_(_linear), constant_(_constant) {}
        //! Performs comparison.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const
            { return math::eq(_pos, reference_, linear_ * (_pos-reference_).norm() + constant_); }
        //! Performs comparison.
        bool operator()(math::rVector3d const& _pos) const
          { return math::eq(_pos, reference_, linear_ * (_pos-reference_).norm() + constant_); }
      protected:
        //! Object to which to compare.
        math::rVector3d reference_;
        //! Linear tolerance, with respect to norm. 
        types::t_real linear_;
        //! constant tolerance
        types::t_real constant_;
    };

    //! Compares both position and types.
    template<class T_TYPE>
      struct CompareSites : public CompareOccupations<T_TYPE>, public ComparePositions
      {
        //! Constructor.
        CompareSites   (Atom<T_TYPE> const _atom, types::t_real _tol=types::tolerance) 
                     : CompareOccupations<T_TYPE>(_atom.type), ComparePositions(_atom.pos, _tol) {}
        template<class T>
        CompareSites   ( Eigen::DenseBase<T> const &_pos,
                         T_TYPE const &_type, 
                         types::t_real _tol=types::tolerance ) 
                     : CompareOccupations<T_TYPE>(_type), ComparePositions(_pos, _tol) {}
        CompareSites(CompareSites const &_c) : CompareOccupations<T_TYPE>(_c), ComparePositions(_c) {}
        //! Compares both position and types.
        bool operator()(Atom<T_TYPE> const& _atom) const
          { return     CompareOccupations<T_TYPE>::operator()(_atom.type)
                   and ComparePositions::operator()(_atom.pos); }
        //! Performs comparison.
        bool operator()( T_TYPE const& _type ) const
          { return CompareOccupations<T_TYPE>::operator()(_type); }
        //! Performs comparison.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const
            { return ComparePositions::operator()(_pos); }
      };

    //! Compares both position and types.
    template<class T_TYPE>
      struct RelCompareSites : public CompareOccupations<T_TYPE>, public RelComparePositions
      {
        //! Constructor.
        RelCompareSites   ( Atom<T_TYPE> const _atom,
                            types::t_real _linear=types::tolerance, 
                            types::t_real _constant=types::tolerance ) 
                        : CompareOccupations<T_TYPE>(_atom.type), 
                          RelComparePositions(_atom.pos, _linear, _constant) {}
        template<class T>
        RelCompareSites   ( Eigen::DenseBase<T> const &_pos,
                            T_TYPE const &_type, 
                            types::t_real _linear=types::tolerance,
                            types::t_real _constant=types::tolerance ) 
                        : CompareOccupations<T_TYPE>(_type), RelComparePositions(_pos, _linear, _constant) {}
        //! Compares both position and types.
        bool operator()(Atom<T_TYPE> const& _atom) const
          { return     CompareOccupations<T_TYPE>::operator()(_atom.type)
                   and RelComparePositions::operator()(_atom.pos); }
        //! Performs comparison.
        bool operator()( T_TYPE const& _type ) const
          { return CompareOccupations<T_TYPE>::operator()(_type); }
        //! Performs comparison.
        template<class T_DERIVED>
          bool operator()(Eigen::DenseBase<T_DERIVED> const& _pos) const
            { return RelComparePositions::operator()(_pos); }
      };


    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T_TYPE>
      CompareSites<T_TYPE> compare_sites( Atom<T_TYPE> const &_origin,
                                          types::t_real _tolerance = types::tolerance )
        { return CompareSites<T_TYPE>(_origin, _tolerance); }
    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T, class T_TYPE>
      CompareSites<T_TYPE> compare_sites( Eigen::DenseBase<T> const &_pos, 
                                          T_TYPE const &_type, 
                                          types::t_real _tolerance = types::tolerance )
        { return CompareSites<T_TYPE>(_pos, _type, _tolerance); }
    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T, class T_TYPE>
      RelCompareSites<T_TYPE> relcompare_sites( Eigen::DenseBase<T> const &_pos, 
                                                T_TYPE const &_type, 
                                                types::t_real _linear = types::tolerance,
                                                types::t_real _constant = types::tolerance )
        { return RelCompareSites<T_TYPE>(_pos, _type, _linear, _constant); }
    //! Returns a functor to compare sites according to positions.
    template<class T_TYPE>
      ComparePositions compare_positions( Eigen::DenseBase<T_TYPE> const &_origin,
                                          types::t_real _tolerance = types::tolerance )
        { return ComparePositions(_origin, _tolerance); }
    //! Returns a functor to compare sites according to positions and/or occupations.
    template<class T_TYPE>
      CompareOccupations<T_TYPE> compare_occupations(T_TYPE const &_origin)
        { return CompareOccupations<T_TYPE>(_origin); }
  } // namespace Crystal
} // namespace LaDa
#endif
