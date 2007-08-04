#ifndef _PESCAN_H_
#define _PESCAN_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string>
#include <string>

#include <eo/eoOp.h>

#include <tinyxml/tinyxml.h>

#include "vff/functional.h"
#include "pescan_interface/interface.h"
#include "lamarck/structure.h"
#include "opt/opt_function_base.h"
#include "opt/opt_minimize_gsl.h"
#include "opt/types.h"

#include "two_sites.h"
#include "evaluator.h"
#include "concentration.h"
#include "functors.h"

#ifdef _MPI
#include "mpi/mpi_object.h"
#endif

namespace BandGap
{
  // actual vff cum pescan functional
  class Functional : public function :: Base<>
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Functional>( BandGap::Functional & );
#endif
    public:
      typedef types::t_real t_Type;
      typedef std::vector<t_Type> t_Container;
      typedef t_Container :: iterator iterator;
      typedef t_Container :: const_iterator const_iterator;

    protected:
      Ising_CE::Structure &structure;
      Vff::Functional vff;
      Pescan::Interface pescan;
      minimizer::GnuSL<Vff::Functional> vff_minimizer;
      std::string references_filename;
      types::t_int nbeval;

    public:
      Functional   ( Ising_CE::Structure &_str, std::string _f = "BandEdge" ) 
                 : structure( _str ), vff(structure),
                   vff_minimizer(vff), references_filename(_f), nbeval(0) {}
      Functional   ( const Functional &_func )
                 : structure(_func.structure), 
                   vff(_func.vff), pescan(_func.pescan), vff_minimizer(_func.vff_minimizer),
                   references_filename(_func.references_filename), nbeval() {};
      ~Functional() {};
      bool Load( const TiXmlElement &_node );
      
      void set_filename( const std::string &_f )
        { references_filename = _f; } 
      void read_references();
      void write_references();
      void set_all_electron() { pescan.set_method( Pescan::Interface::Escan::ALL_ELECTRON ); }
      t_Type evaluate()
      {
        // first minimizes strain
        std::ostringstream sstr; sstr << "escan" << nbeval; 
        ++nbeval;
#ifdef _MPI
        sstr << mpi::main.rank();
#endif
        pescan.set_dirname(sstr.str());
        vff_minimizer.minimize();
        structure.energy = vff.energy();
        vff.print_escan_input();

        // then evaluates band gap
        types::t_real result = pescan(structure);
        if ( pescan.get_method() != Pescan::Interface::Escan::ALL_ELECTRON )
          return result;
        
        pescan.set_method(); // resets to folded spectrum if necessary

#ifdef _MPI
        if ( mpi::main.rank() != mpi::ROOT_NODE ) // not root no read write
          return result;
#endif 
        if ( result < types::tolerance ) 
        {
          set_all_electron();
          evaluate();
        }
        write_references();

        return result;
      } 
      void evaluate_gradient( t_Type* const _i_grad ) 
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
      t_Type evaluate_one_gradient( types::t_unsigned _pos)
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
      t_Type evaluate_with_gradient( t_Type* const _i_grad )
      { 
        std::cerr << "Not implemented !!" << std::endl;
        exit(0);
      }
      void get_bands( types::t_real &_vbm, types::t_real &_cbm ) 
        { pescan.get_bands( _vbm, _cbm); }
  };

  struct Object : public TwoSites::Object
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Object>(BandGap::Object &);
#endif
    types::t_real CBM, VBM;

    Object() : TwoSites::Object(), CBM(0), VBM(0) {}
    Object(const Object &_c) : TwoSites::Object(_c), CBM(_c.CBM), VBM(_c.VBM) {};
    ~Object() {};
    
  };


  class Evaluator : public TwoSites::Evaluator< individual::Single<TwoSites::Object> >
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<BandGap::Evaluator>( BandGap::Evaluator & );
#endif
    protected:
      typedef TwoSites::Evaluator< individual::Single<TwoSites::Object> > t_Base;
      typedef Ising_CE::Structure::t_kAtoms t_kvecs;
      typedef Ising_CE::Structure::t_Atoms t_rvecs;

    protected:
      using t_Base :: Load;
      
    protected:
      Vff::Functional vff;
      Pescan::Interface pescan;
      minimizer::GnuSL<Vff::Functional> vff_minimizer;
      std::string references_filename;
      types::t_int nbeval;
      types::t_int age, check_ref_every;

    public:
      Evaluator() : TwoSites::Evaluator<Object>(), 
                    vff(structure), vff_minimizer(vff), references_filename(_f), 
                    nbeval(0)age(0), check_ref_every(-1) {}
      ~Evaluator() {};

      bool initialize( t_Individual &_indiv )
        { return t_Base::initialize( _indiv ); }
      void init( t_Individual &_indiv );
      bool Load( const TiXmlElement &_node );
      bool Load ( t_Individual &_indiv, const TiXmlElement &_node, bool _type );
      bool Save ( const t_Individual &_indiv, TiXmlElement &_node, bool _type ) const;
      void LoadAttribute ( const TiXmlAttribute &_att ) {};
      eoF<bool>* LoadContinue(const TiXmlElement &_el )
        { return new darwin::mem_zerop_t<Evaluator>( *this, &Evaluator::Continue,
                                                     "Evaluator::Continue" );     }

      bool Continue();
      void evaluate( )

    protected:
      void set_all_electron() { pescan.set_method( Pescan::Interface::Escan::ALL_ELECTRON ); }
      void read_references();
      void write_references();
      void get_bands( types::t_real &_vbm, types::t_real &_cbm ) 
        { pescan.get_bands( _vbm, _cbm); }
      void set_filename( const std::string &_f )
        { references_filename = _f; } 
  };

  void rearrange_structure(Ising_CE::Structure &);


} // namespace BandGap

#endif // _PESCAN_H_
