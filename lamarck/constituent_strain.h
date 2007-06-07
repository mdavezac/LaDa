#ifndef _CONSTITTUENT_STRAIN_H_
#define _CONSTITTUENT_STRAIN_H_

#include <vector>
#include <tinyxml/tinyxml.h>
#include <opt/opt_function_base.h>
#include <opt/types.h>

#include "atat/vectmac.h"

#include "structure.h"
#include "harmonic.h"

#ifdef _MPI 
  #include<mpi/mpi_object.h>
#endif

namespace Ising_CE 
{

  class Constituent_Strain : public function::Base<types::t_real>
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Constituent_Strain> ( Constituent_Strain& );
#endif
    public: 
      typedef std::vector<Harmonic> t_Harmonics;

    protected:
      static const types::t_real ZERO_TOLERANCE;

    protected: 
      std::vector<atat::rVector3d> r_vecs, k_vecs;
      atat::rMatrix3d r_cell, k_cell;
      static t_Harmonics harmonics;
 
      // constructor, destructor, and helpers
    public:
      Constituent_Strain() : function::Base<types::t_real>() {};
      Constituent_Strain(const Ising_CE::Structure& str, const atat::rMatrix3d &lattice)
        : function::Base<types::t_real>()
          { set_structure( str, lattice ); }
      Constituent_Strain(const Ising_CE::Structure& str, const atat::rMatrix3d &lattice,
                         std::vector<types::t_real> *vars)
        { function::Base<types::t_real>::variables = vars; set_structure( str, lattice ); }
      ~Constituent_Strain(){};
      void set_structure( const Ising_CE::Structure& str, const atat::rMatrix3d &lattice );
   
      // required behaviors for interfacing with minimizer
    public: 
      types::t_real evaluate();
      void evaluate_gradient(types::t_real* const gradient)
        { evaluate_with_gradient( gradient ); }
      types::t_real evaluate_with_gradient(types::t_real* const gradient);
      types::t_real evaluate_one_gradient( types::t_unsigned _pos );

      // others
    protected:
      void cut_integer_part( atat::rVector3d &kvec);
      void find_range( const atat::rMatrix3d &A, atat::iVector3d &kvec);
      void refold( atat::rVector3d &vec, const atat::rMatrix3d &lat );
      bool are_equivalent( const atat::rVector3d &vec_a,
                           const atat::rVector3d &vec_b,
                           const atat::rMatrix3d &lat) const;
      void remove_equivalents( const atat::rMatrix3d &lat);



    public:
      bool Load_Harmonics( const TiXmlElement &_element);
      bool Load (const TiXmlElement &_element);
      void print_xml( TiXmlElement& _node ) const;
      const std::vector<atat::rVector3d>& get_kvectors() const
          { return k_vecs; }

      #ifdef _DEBUG_LADA_
        void check_derivative();
      #endif // _DEBUG_LADA_

  };

} // namespace Ising_CE
#endif // _CONSTITTUENT_STRAIN_H_
