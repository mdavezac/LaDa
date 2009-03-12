//
//  Version: $Id$
//
#ifndef _PESCAN_EMASS_H_
#define _PESCAN_EMASS_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <utility>

#include <minimizer/cgs.h>
#include <crystal/structure.h>
#include <atat/vectmac.h>
#include <opt/types.h>

#include "interface.h"

namespace LaDa
{
  namespace Pescan 
  {
    class eMass
    {
      public:
        //! Conjugate gradient minimizer.
        typedef Fitting :: Cgs t_Minimizer;
        //! Conjugate gradient minimizer.
        t_Minimizer cgs;
        //! Order of interpolation.
        size_t order;
        //! Number of interpolation points.
        size_t npoints;
        //! Stepsize for interpolation points.
        types::t_real stepsize;

        //! Constructor.
        eMass() : order( 2 ), npoints(2), stepsize( 1e-2 ) {}
        //! Copy Constructor.
        eMass   ( const eMass &_c )
              : cgs( _c.cgs ), order( _c.order ), 
                npoints( _c.npoints ), stepsize( _c.stepsize ) {}

        //! \brief Computes effective mass.
        //! \param[in] _interface the escan functional.
        //! \param[in] _ocell cell vectors of the structure prior to relaxation.
        //! \param[in] _structure \e relaxed structure for which to compute emasses.
        //! \param[in] _at effective masses at vector \a _at.
        //! \param[in] _direction of effective mass.
        //! \param[in] _nbstates number of bands (and effective masses) to compute.
        //! \param[in] _eref reference energy at which to compute bands.
        //! \param[out] _out resulting effective masses and band eigenvalues.
        void operator()
        (
          const Interface& _interface,
          const atat::rMatrix3d &_ocell, 
          const Crystal::Structure &_structure,
          const atat::rVector3d &_at,
          const atat::rVector3d &_direction,
          const size_t &_nbstates,
          const types::t_real &_eref,
          std::vector< std::pair<types::t_real, types::t_real> > &_out 
        ) const;

      protected:
        //! Computes kpoint in distorted structure from original lattice point.
        atat::rVector3d distort_kpoint( const atat::rMatrix3d &_ocell, 
                                        const atat::rMatrix3d &_dcell,
                                        const atat::rVector3d &_kpoint ) const
          { return ( !(~_dcell) ) * (~_ocell) * _kpoint; }
        //! Computes eigenvalues for single kpoint.
        void compute_( Interface& _interface,
                       const atat::rVector3d &_kpoints ) const;

        //! Adds a weight to interpolation points.
        types::t_real weight_( size_t _i, size_t j ) const;

    };

  } // namespace Pescan
} // namespace LaDa

#endif // _PESCAN_EMASS_H_
