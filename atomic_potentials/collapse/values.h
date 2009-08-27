//
//  Version: $Id: representations.h 1271 2009-08-17 17:57:45Z davezac $
//
#ifndef LADA_ATOMIC_POTENTIAL_COLLAPSE_VALUES_H_
#define LADA_ATOMIC_POTENTIAL_COLLAPSE_VALUES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include <opt/debug.h>
#include <opt/fuzzy.h>

#include "../numeric_types.h"

namespace LaDa
{
# ifdef __DOPYTHON
  //! \cond
  namespace Python
  {
    void expose_values();
  }
  //! \endcond
# endif

  namespace atomic_potential
  {
    // Forward declaration
    //! \cond
    class Representation;
    class SumOfSeparables;
    //! \endcond

    namespace collapse
    {

      // Forward declaration
      //! \cond
      class FittingSet;
      //! \endcond

      //! Function values and aggregates thereof required by alternating least-square fit.
      class Values
      {
        friend void ::LaDa::Python::expose_values();
        protected:
          //! Type of the Fitting set.
          typedef FittingSet t_FittingSet;
          //! \brief Type of the containers of factors (i!=j) for alternating least square fit.
          //! \details elements: \f$crv =  \sum_o \alpha_{o,r,i}[x_i^{(u,v)}] *
          //!                       g_{o,r,i}[x_i^{(u,v)}]\f$.
          typedef std::vector // over coordinates.
                  <
                    std::vector // over structures.
                    <
                      std::vector // over representations.
                      <
                        std::vector<numeric_type> // over ranks.
                      >
                    >
                  > t_CoordRankValues;
          //! \brief Type of the containers of factors (i!=j) for alternating least square fit.
          //! \details elements: \f$rv =  \Prod_i\sum_o \alpha_{o,r,i}[x_i^{(u,v)}] *
          //!             g_{o,r,i}[x_i^{(u,v)}]\f$.
          typedef std::vector // over structures.
                  <
                    std::vector // over representations.
                    <
                      std::vector<numeric_type> // over ranks.
                    >
                  > t_RankValues;
          //! \brief Type of the containers of function values.
          //! \details elements: \$[ g_{o,r,i}[x_i^{(u,v)}]\f$
          typedef std::vector // over coordinates
                  < 
                    std::vector // over structures
                    <
                      std::vector // over representations
                      <
                        std::vector // over ranks
                        <
                          std::vector<numeric_type> // inner sum.
                        >
                      >
                    >
                  > t_FunctionValues;
        public:
          //! An iterator over the structures.
          class str_iterator;
          //! An iterator over the structures.
          class const_str_iterator;
#         include "values.str_iterator.h"
#         ifdef LADA_WITH_CONST
#           error LADA_WITH_CONST already defined.
#         endif
#         define LADA_WITH_CONST
#         include "values.str_iterator.h"
          //! Constructor.
          Values() {}
          //! Copy Constructor.
          Values   (Values const &_c)
                 : coord_rank_values_(_c.coord_rank_values_),
                   rank_values_(_c.rank_values_),
                   function_values_(_c.function_values_) {}
          //! returns iterator to start of structures, for variable \a _i.
          str_iterator begin( size_t _i );
          //! returns iterator to end of structures, for variable \a _i.
          str_iterator end( size_t _i );
          //! returns iterator to start of structures, for variable \a _i.
          const_str_iterator begin( size_t _i ) const;
          //! returns iterator to end of structures, for variable \a _i.
          const_str_iterator end( size_t _i ) const;

          //! Updates values using new coefficients.
          void update( vector_type const& _coefs, t_FittingSet const &_ftstr, size_t _i );
          //! Adds new structure to fitting set.
          void add( Representation const &_rep, SumOfSeparables const &_sep );

        protected:
          //! Coordinates for each structure and representation.
          t_CoordRankValues coord_rank_values_;
          //! Coordinates for each structure and representation.
          t_RankValues rank_values_;
          //! Weight of each structure and represenation.
          t_FunctionValues function_values_;
      };

//     class Values::str_iterator
//     {
//         friend class Values;
//       public:
//         //! Iterator over representations.
//         class rep_iterator;
//         //! Constructor.
//         str_iterator   (size_t _i, t_CoordRankValues const &_cr)
//                      : i_(_i), n_(0), coord_rank_values_(_cr) {}
//         //! Copy Constructor.
//         str_iterator   (str_iterator const &_c) 
//                      : i_(_c.i_), n_(_c.n_), coord_rank_values_(_c.coord_rank_values_),
//                        i_coord_rank_value_(_c.i_coord_rank_value_),
//                        i_rank_value_(_c.i_rank_value_),
//                        i_function_values_(_c.i_function_values_) {}
//
//         //! Increments iterator. 
//         void operator++() 
//         {
//           ++i_coord_rank_value_;
//           ++i_rank_value_;
//           ++i_function_values_;
//           ++n_;
//         }
//
//         //! Decrements iterator.
//         void operator--() 
//         {
//           --i_coord_rank_value_;
//           --i_rank_value_;
//           --i_function_values_;
//           --n_;
//         }
//         //! Returns iterator over representations.
//         rep_iterator begin() const;
//         //! Returns iterator over representations.
//         rep_iterator end() const;
//
//         //  True if iterators are at same position.
//         bool operator==(str_iterator const& _b) const
//         {
//           LADA_ASSERT( i_ == _b.i_, "Inequivalent iterators.\n");
//           return i_coord_rank_value_ != _b.i_coord_rank_value_; 
//         }
//         //! False if iterators are at same position.
//         bool operator!=(str_iterator const& _b) const { return not operator==(_b); }
//
//       protected:
//         //! Current variable.
//         size_t i_;
//         //! Current variable.
//         size_t n_;
//         //! Values aggregate per coordinate and rank.
//         t_CoordRankValues const& coord_rank_values_;
//         //! Iterator over other factors.
//         t_CoordRankValues::value_type::const_iterator i_coord_rank_value_;
//         //! Iterator over other factors.
//         t_RankValues::const_iterator i_rank_value_;
//         //! Iterator over function values.
//         t_FunctionValues::value_type::const_iterator i_function_values_;
//     };
//
//
//     class Values::str_iterator::rep_iterator
//     {
//         friend class Values::str_iterator;
//       public:
//         //! Iterator over ranks.
//         class rank_iterator; 
//         //! Constructor.
//         rep_iterator   (size_t _i, size_t _nstr, t_CoordRankValues const &_cr)
//                      : i_(_i), n_str_(_nstr), n_(0), coord_rank_values_(_cr) {}
//         //! Copy Constructor.
//         rep_iterator   (rep_iterator const &_c) 
//                      : i_(_c.i_), n_str_(_c.n_str_), n_(_c.n_),
//                        coord_rank_values_(_c.coord_rank_values_),
//                        i_coord_rank_value_(_c.i_coord_rank_value_),
//                        i_rank_value_(_c.i_rank_value_),
//                        i_function_values_(_c.i_function_values_) {}
//
//         //! Returns iterator over representations.
//         rank_iterator begin() const;
//         //! Returns iterator over representations.
//         rank_iterator end() const;
//
//         //! Increments operator.
//         void operator++() { ++i_coord_rank_value_; ++i_rank_value_; ++i_function_values_; ++n_; }
//         //! Decrements operator.
//         void operator--() { --i_coord_rank_value_; --i_rank_value_; --i_function_values_; --n_; }
//         //  True if iterators are at same position.
//         bool operator==(rep_iterator const& _b) const
//         {
//           LADA_ASSERT( i_ == _b.i_, "Inequivalent iterators.\n");
//           return i_coord_rank_value_ != _b.i_coord_rank_value_; 
//         }
//         //! False if iterators are at same position.
//         bool operator!=(rep_iterator const& _b) const { return not operator==(_b); }
//         
//       protected:
//         //! Current variable.
//         size_t i_;
//         //! Current structure.
//         size_t n_str_;
//         //! Current representation.
//         size_t n_;
//         //! Values aggregate per coordinate and rank.
//         t_CoordRankValues const& coord_rank_values_;
//         //! Iterator over other factors.
//         t_CoordRankValues::value_type::value_type::const_iterator i_coord_rank_value_;
//         //! Iterator over other factors.
//         t_RankValues::value_type::const_iterator i_rank_value_;
//         //! Iterator over function values.
//         t_FunctionValues::value_type::value_type::const_iterator i_function_values_;
//     };
//
//
//
//     class Values::str_iterator::rep_iterator::rank_iterator
//     {
//         friend class Values::str_iterator::rep_iterator;
//       public:
//         //! Iterator over function values.
//         typedef t_FunctionValues::value_type::value_type::value_type
//                                 ::value_type::const_iterator function_iterator;
//         //! Constructor.
//         rank_iterator  (size_t _i, size_t _nstr, size_t _nrep, t_CoordRankValues const& _cr)
//                       : i_(_i), n_str_(_nstr), n_rep_(_nrep), n_(0), coord_rank_values_(_cr) {}
//         //! Copy Constructor.
//         rank_iterator   (rank_iterator const &_c) 
//                       : i_(_c.i_), n_str_(_c.n_str_), n_rep_(_c.n_rep_), n_(_c.n_),
//                         coord_rank_values_(_c.coord_rank_values_),
//                         i_coord_rank_value_(_c.i_coord_rank_value_),
//                         i_rank_value_(_c.i_rank_value_),
//                         i_function_values_(_c.i_function_values_) {}
//
//         //! Returns factors (i!=j) for alternating least-square-fit.
//         numeric_type other() const;
//         //! Returns iterator to function values.
//         function_iterator begin() const { return i_function_values_->begin(); }
//         //! Returns iterator to function values.
//         function_iterator end() const { return i_function_values_->end(); }
//         //! Increments operator.
//         void operator++() 
//           { ++i_coord_rank_value_; ++i_rank_value_; ++i_function_values_; ++n_; }
//
//
//        //  True if iterators are at same position.
//        bool operator==(rank_iterator const& _b ) const
//          { return i_coord_rank_value_ != _b.i_coord_rank_value_;  }
//        //! False if iterators are at same position.
//        bool operator!=(rank_iterator const& _b ) const { return not operator==(_b); }
//       protected:
//         //! Current variable.
//         size_t i_;
//         //! Current structure.
//         size_t n_str_;
//         //! Current representation.
//         size_t n_rep_;
//         //! Current rank
//         size_t n_;
//         //! Values aggregate per coordinate and rank.
//         t_CoordRankValues const& coord_rank_values_;
//         //! Iterator over other factors.
//         t_CoordRankValues::value_type::value_type::value_type::const_iterator i_coord_rank_value_;
//         //! Iterator over other factors.
//         t_RankValues::value_type::value_type::const_iterator i_rank_value_;
//         //! Iterator over function values.
//         t_FunctionValues::value_type::value_type::value_type::const_iterator i_function_values_;
//     };

    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif

