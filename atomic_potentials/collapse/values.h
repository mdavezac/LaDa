//
//  Version: $Id: representations.h 1271 2009-08-17 17:57:45Z davezac $
//
#ifndef LADA_ATOMIC_POTENTIAL_COLLAPSE_VALUES_H_
#define LADA_ATOMIC_POTENTIAL_COLLAPSE_VALUES_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>

#include "../numeric_types.h"

namespace LaDa
{
  namespace atomic_potential
  {
    namespace collapse
    {
      class FittingSet;

      //! Function values and aggregates thereof required by alternating least-square fit.
      class Values
      {
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
                  > t_CoordRankValue;
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
          //! Constructor.
          Values() {}
          //! Copy Constructor.
          Values   (Values const &_const)
                 : coord_rank_values_(_c.coord_rank_values_),
                   rank_values_(_c.rank_values_),
                   function_values_(_c.function_values_) {}
          //! returns iterator to start of structures, for variable \a _i.
          str_iterator begin( size_t _i ) const;
          //! returns iterator to end of structures, for variable \a _i.
          str_iterator end( size_t _i ) const;

          //! Updates values using new coefficients.
          void update( vector_type const& _coefs, t_FittingSet const &_ftstr, size_t _i );
          //! Adds new structure to fitting set.
          void add(t_FittingSet const &_ftstr, Crystal::TStructure<std::string> _structure);

        protected:
          //! Coordinates for each structure and representation.
          t_CoordRankValue coord_rank_values_;
          //! Coordinates for each structure and representation.
          t_RankValue coord_values_;
          //! Weight of each structure and represenation.
          t_FunctionValues function_values_;
      };

      //! True if iterators are at same position.
      bool operator==( Values::str_iterator const& _a,
                       Values::str_iterator const& _b );

      class Values::str_iterator
      {
          friend bool operator==( str_iterator const& _a, str_iterator const& _b );
        public:
          //! Iterator over representations.
          class rep_iterator;
          //! Constructor.
          str_iterator   (size_t _i, t_CoordRankValues const &_cr)
                       : i_(_i), n_(0), coord_rank_values_(_cr) {}
          //! Copy Constructor.
          str_iterator   (str_iterator const &_c) 
                       : i_(_c.i_), n_(_c.n_), coord_rank_values_(_c.coord_rank_values_),
                         i_coord_rank_value_(_c.i_coord_rank_value_),
                         i_rank_value_(_c.i_rank_value_),
                         i_function_values_(_c.i_function_values_) {}

          //! Increments iterator. Returns true if not end of container.
          void operator++() 
          {
            ++i_coord_rank_value_;
            ++i_rank_value_;
            ++i_function_values_;
          }

          //! Returns iterator over representations.
          rep_iterator begin() const;
          //! Returns iterator over representations.
          rep_iterator end() const;

        protected:
          //! Current variable.
          size_t i_;
          //! Current variable.
          size_t n_;
          //! Values aggregate per coordinate and rank.
          t_CoordRankValue const& coord_rank_values_;
          //! Iterator over other factors.
          t_CoordRankValue::value_type::const_iterator i_coord_rank_value_;
          //! Iterator over other factors.
          t_RankValue::const_iterator i_rank_value_;
          //! Iterator over function values.
          t_FunctionValues::value_type::const_iterator i_function_values_;
      };

      //  True if iterators are at same position.
      inline bool operator==( Values::str_iterator const& _a,
                              Values::str_iterator const& _b )
      {
        LADA_ASSERT( _a.i_ == _b.i_, "Inequivalent iterators.\n");
        return _a.i_coord_rank_value_ != _b.i_coord_rank_value_; 
      }
      //! False if iterators are at same position.
      inline bool operator!=( Values::str_iterator const& _a,
                              Values::str_iterator const& _b )
        { return not (_a == _b); }

      //! True if iterators are at same position.
      bool operator==( Values::str_iterator::rep_iterator const& _a,
                       Values::str_iterator::rep_iterator const& _b );
      
      class Values::str_iterator::rep_iterator
      {
          friend bool operator==( rep_iterator const& _a, rep_iterator const& _b );
        public:
          //! Iterator over ranks.
          class rank_iterator; 
          //! Constructor.
          rep_iterator   (size_t _i, size_t _nstr, t_CoordRankValues const &_cr)
                       : i_(_i), n_str_(_nstr), n_(0), coord_rank_values_(_coord_rank_values_) {}
          //! Copy Constructor.
          rep_iterator   (rep_iterator const &_c) 
                       : i_(_c.i), n_str_(_c.n_str_), n_(_c.n_),
                         coord_rank_values_(_c.coord_rank_values_),
                         i_coord_rank_value_(_c.i_coord_rank_value_),
                         i_rank_value_(_c.i_rank_value_),
                         i_function_values_(_c.i_function_values_) {}

          //! Returns iterator over representations.
          rank_iterator begin() const;
          //! Returns iterator over representations.
          rank_iterator end() const;
          
        protected:
          //! Current variable.
          size_t i_;
          //! Current structure.
          size_t n_str_;
          //! Current representation.
          size_t n_;
          //! Values aggregate per coordinate and rank.
          t_CoordRankValue const& coord_rank_values_;
          //! constant reference to container of coefficients.
          t_Coefficients &coefficients_;
          //! Iterator over other factors.
          t_CoordRankValue::value_type::value_type::const_iterator i_coord_rank_value_;
          //! Iterator over other factors.
          t_RankValue::value_type::const_iterator i_rank_value_;
          //! Iterator over function values.
          t_FunctionValues::value_type::value_type::const_iterator i_function_values_;
      };

      //  True if iterators are at same position.
      inline bool operator==( Values::str_iterator::rep_iterator const& _a,
                              Values::str_iterator::rep_iterator const& _b )
      {
        LADA_ASSERT( _a.i_ == _b.i_, "Inequivalent iterators.\n");
        return _a.i_coord_rank_value_ != _b.i_coord_rank_value_; 
      }
      //! False if iterators are at same position.
      inline bool operator!=( Values::str_iterator::rep_iterator const& _a,
                              Values::str_iterator::rep_iterator const& _b )
        { return not (_a == _b); }

      //! True if iterators are at same position.
      bool operator==( Values::str_iterator::rep_iterator::rank_iterator const& _a,
                       Values::str_iterator::rep_iterator::rank_iterator const& _b );

      class Values::str_iterator::rep_iterator::rank_iterator
      {
        friend bool operator==( rank_iterator const& _a, rank_iterator const& _b );
        public:
          //! Iterator over function values.
          typedef t_FunctionValues::value_type::value_type
                                  ::value_type::const_iterator function_iterator;
          //! Constructor.
          rank_iterator  (size_t _i, size_t _nstr, size_t _nrep, t_CoordRankValues const& _cr)
                        : i_(_i), n_str_(_nstr), n_rep_(_nrep), n_(0), coord_rank_values_(_cr) {}
          //! Copy Constructor.
          rank_iterator   (rank_iterator const &_c) 
                        : i_(_c.i_), n_str_(_c.n_str_), n_rep_(_c.n_rep_), n_(_c.n_),
                          coord_rank_values_(_c.coord_rank_values_),
                          i_coord_rank_value_(_c.i_coord_rank_value_),
                          i_rank_value_(_c.i_rank_value_),
                          i_function_values_(_c.i_function_values_) {}

          //! Returns factors (i!=j) for alternating least-square-fit.
          numeric_type other() const;
          //! Returns iterator to function values.
          function_iterator begin() const { return i_function_values_->begin(); }
          //! Returns iterator to function values.
          function_iterator end() const { return i_function_values_->end(); }

        protected:
          //! Current variable.
          size_t i_;
          //! Current structure.
          size_t n_str_;
          //! Current representation.
          size_t n_rep_;
          //! Current rank
          size_t n_;
          //! Values aggregate per coordinate and rank.
          t_CoordRankValue const& coord_rank_values_;
          //! Iterator over other factors.
          t_CoordRankValue::value_type::value_type::value_type::const_iterator i_coord_rank_value_;
          //! Iterator over other factors.
          t_RankValue::value_type::value_type::const_iterator i_rank_value_;
          //! Iterator over function values.
          t_FunctionValues::value_type::value_type::const_iterator i_function_values_;
      };

      //  True if iterators are at same position.
      inline bool operator==( Values::str_iterator::rep_iterator::rank_iterator const& _a,
                              Values::str_iterator::rep_iterator::rank_iterator const& _b )
        { return _a.i_coord_rank_value_ != _b.i_coord_rank_value_;  }
      //! False if iterators are at same position.
      inline bool operator!=( Values::str_iterator::rep_iterator::rank_iterator const& _a,
                              Values::str_iterator::rep_iterator::rank_iterator const& _b )
        { return not (_a == _b); }

    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif

