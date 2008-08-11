//
//  Version: $Id$
//
#ifndef _CE_COLLAPSE_H_
#define _CE_COLLAPSE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include<boost/numeric/ublas/vector.hpp>
#include<boost/numeric/ublas/matrix.hpp>

#include <opt/types.h>
#include <opt/debug.h>

#include "prepare.h"

namespace CE
{
  namespace Methods
  {
    template< class T_COLLAPSE, class T_SEPARABLES, class T_STRUCTURES >
      opt::ErrorTuple check_one( const T_SEPARABLES &_separables,
                                 const T_COLLAPSE &_collapse,
                                 const T_STRUCTURES &_str,
                                 size_t _n, bool _verbose = false );
  }
                                                   
  //! Collapse functor for fitting CE::Separables  
  template< class T_SEPARABLES, class T_MAPPING,
            class T_NORMALIZATION = typename T_SEPARABLES :: t_Mapping >
    class Collapse
    {
      template< class T_COLLAPSE, class TT_SEPARABLES, class T_STRUCTURES > friend
        opt::ErrorTuple Methods::check_one( const TT_SEPARABLES &_separables,
                                            const T_COLLAPSE &_collapse,
                                            const T_STRUCTURES &_str,
                                            size_t _n, bool _verbose );
      public:
        //! Type of the separable function.
        typedef T_SEPARABLES t_Separables;
        //! Type of the mapping function from structures to targets.
        typedef T_MAPPING t_Mapping;
        //! Type of normalization used.
        typedef T_NORMALIZATION t_Normalization;
        //! Type of the matrices.
        typedef typename t_Separables :: t_Matrix t_Matrix;
        //! Type of the vectors.
        typedef typename t_Separables :: t_Vector t_Vector;

        //! \brief The mapping from target values to symetrically equivalent
        //!        structures.
        t_Mapping mapping;

        //! Constructor.
        Collapse() : dim(0), separables(NULL) {}
        //! Destructor.
        ~Collapse() {}

        //! Creates the fitting matrix and target vector.
        void operator()( t_Matrix &_A, t_Vector &_b,
                         types::t_unsigned &_dim )
          { dim = _dim; create_A_n_b( _A, _b ); }
        //! Evaluates square errors.
        opt::ErrorTuple evaluate();

        //! Updates the scales vector and  normalizes.
        void update_all();
        //! Updates the scales vector, should update only one dim.
        void update( types::t_unsigned _d )
          { update_all(); }
        //! Does nothing.
        void reset() {}

        //! Initializes collapse functor.
        template< class T_STRUCTURES >
        void init( const T_STRUCTURES& _strs, const PosToConfs &_postoconfs );
        //! Sets the norm pointer.
        void init( t_Separables& _sep ) { separables = &_sep; }

        //! Reference to configuration matrix.
        const t_Matrix& configurations() const { return configurations_; }

      protected:
        //! Creates the _A and _b matrices for fitting.
        void create_A_n_b( t_Matrix &_A, t_Vector &_b );
        //! Finds scaling factor for that conf, collapsed dimension, and rank.
        typename t_Vector::value_type factor( size_t _kv, size_t _r );
        //! \brief Creates an X vector for fitting, for a single rank.
        //! \details Before the scaling is done, all ranks equal.
        void create_X( size_t _i, t_Vector &_out );


        //! The configurations, arranged in columns.
        t_Matrix configurations_;
        //! Holds current dimension being fitted.
        size_t dim;
        //! holds the sep. func. split up by rank and confs.
        t_Matrix scales;
        //! Pointer to separable function being minimized.
        t_Separables *separables;
    };

}

#include "collapse.impl.h"

#endif
