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
          void update(vector_type const& _coefs, t_FittingSet const &_ftstr, size_t _i );
          //! Adds new structure to fitting set.
          void add( Representation const &_rep, SumOfSeparables const &_sep );

          //! Number of coordinates.
          size_t nb_coordinates() const { return coord_rank_values_.size(); }

        protected:
          //! Coordinates for each structure and representation.
          t_CoordRankValues coord_rank_values_;
          //! Coordinates for each structure and representation.
          t_RankValues rank_values_;
          //! Weight of each structure and represenation.
          t_FunctionValues function_values_;
      };


    } // namespace collapse
  } // namespace atomic_potential
} // namespace LaDa
#endif

