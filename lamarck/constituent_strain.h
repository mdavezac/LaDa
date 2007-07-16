#ifndef _CONSTITTUENT_STRAIN_H_
#define _CONSTITTUENT_STRAIN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <tinyxml/tinyxml.h>
#include <opt/opt_function_base.h>
#include <opt/types.h>

#include "atat/vectmac.h"

#include "structure.h"
#include "harmonic.h"

#ifdef _MPI 
  #include "mpi/mpi_object.h"
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
      typedef types::t_real t_Type;
      typedef std::vector<types::t_real>  t_Container;

    protected:
      static const types::t_real ZERO_TOLERANCE;

    protected: 
      std::vector<atat::rVector3d> r_vecs, k_vecs;
      static t_Harmonics harmonics;
 
      // constructor, destructor, and helpers
    public:
      Constituent_Strain() : function::Base<types::t_real>() {};
      Constituent_Strain(const Ising_CE::Structure& str, t_Container *vars=NULL);
      ~Constituent_Strain(){};

   
      // required behaviors for interfacing with minimizer
    public: 
      types::t_real evaluate();
      void evaluate_gradient(types::t_real* const gradient)
        { evaluate_with_gradient( gradient ); }
      types::t_real evaluate_with_gradient(types::t_real* const gradient);
      types::t_real evaluate_one_gradient( types::t_unsigned _pos );

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
