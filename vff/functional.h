#ifndef _VFF_FUNCTIONAL_H_
#define _VFF_FUNCTIONAL_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <vector>
#include <algorithm>

#include <tinyxml/tinyxml.h>

#include "lamarck/atom.h"
#include "lamarck/structure.h"
#include "lamarck/lattice.h"
#include "opt/types.h"
#include "opt/opt_function_base.h"

#ifdef _MPI 
  #include "mpi/mpi_object.h"
#endif

namespace Vff
{

  class Functional;

  class Atomic_Center
  {
    friend class Functional;
    public:
      class const_iterator;

    protected:
      Ising_CE::Atom *origin;
      std::vector< Atomic_Center* > bonds;
      std::vector< atat::rVector3d > translations;
      std::vector< bool > do_translates;
      Ising_CE :: Structure *structure;
      bool is_site_one, is_site_one_two_species;
      atat::rVector3d gradient;
      types::t_unsigned index;

    public:
      Atomic_Center ( Ising_CE::Structure &_str, Ising_CE::Atom &_e, types::t_unsigned _i);
      Atomic_Center   ( const Atomic_Center &_c )
                    : origin(_c.origin), bonds(_c.bonds), translations(_c.translations), 
                      do_translates(_c.do_translates), structure(_c.structure),
                      is_site_one(_c.is_site_one),
                      is_site_one_two_species( _c.is_site_one_two_species), 
                      gradient(0,0,0), index( _c.index) {} 

      types::t_unsigned kind() const;
      types::t_int add_bond( Atomic_Center &_e, const types::t_real _cutoff  );

      const_iterator begin() const;
      const_iterator end() const;
      types::t_unsigned size() const
        { return bonds.size(); }

      atat::rVector3d& operator=(atat::rVector3d& _vec)
        { origin->pos = _vec; return _vec; }
      void operator+=(const atat::rVector3d& _vec)
        { origin->pos += _vec; }
      void operator-=(const atat::rVector3d& _vec)
        { origin->pos -= _vec; }
      operator atat::rVector3d& ()
        { return origin->pos; }
      operator const atat::rVector3d& () const
        { return origin->pos; }
      Ising_CE::Atom& get_origin()
        { return *origin; }
      const Ising_CE::Atom& get_origin() const
        { return *origin; }
      atat::rVector3d& get_gradient()
        { return gradient; }
      const atat::rVector3d& get_gradient() const
        { return gradient; }
      void reset_gradient()
        { gradient[0] = 0; gradient[1] = 0; gradient[2] = 0; }
      bool site_one() const
        { return is_site_one; }
      types::t_unsigned get_index() const
        { return index; }

    protected:
      types::t_unsigned bond_kind( const Atomic_Center &_bond ) const;
  };

  class Atomic_Center :: const_iterator
  {
    protected:
      const Atomic_Center *parent;
      std::vector< Atomic_Center* > :: const_iterator i_bond;
      std::vector< atat::rVector3d > :: const_iterator i_translation;
      std::vector< bool > :: const_iterator i_do_translate;
    
    public:
      const_iterator() {};
      const_iterator   ( const Atomic_Center *_parent, bool _is_begin=true) 
                     : parent(_parent)
      {
        if ( _is_begin )
        {
          i_bond = parent->bonds.begin();
          i_translation = parent->translations.begin();
          i_do_translate = parent->do_translates.begin();
          return;
        }
        i_bond = parent->bonds.end();
      }
      const_iterator   ( const const_iterator &_c ) 
                     : parent(_c.parent), i_bond(_c.i_bond), 
                       i_translation(_c.i_translation),
                       i_do_translate(_c.i_do_translate) {};
      void operator ++()
        { ++i_bond; ++i_translation; ++i_do_translate; }
      void operator --()
        { --i_bond; --i_translation; --i_do_translate; }
      void operator -=( types::t_int  n)
        { i_bond -= n; i_translation -= n; i_do_translate -= n; }
      void operator +=( types::t_int  n)
        { i_bond += n; i_translation += n; i_do_translate += n; }
      bool operator ==( const const_iterator &_i ) const
        { return _i.i_bond == i_bond; }
      bool operator !=( const const_iterator &_i ) const
        { return _i.i_bond != i_bond; }
      types::t_int operator -( const const_iterator &_i )
        { return _i.i_bond - i_bond; }
      Atomic_Center& operator *()
        { return *(*i_bond); }
      const Atomic_Center& operator *() const
        { return *(*i_bond); }
      Atomic_Center* operator ->()
        { return *i_bond; }
      const Atomic_Center* operator ->() const
        { return *i_bond; }
      types::t_real norm2() const
      {
        if ( not *i_do_translate )
          return atat::norm2 ( parent->origin->pos - (*i_bond)->origin->pos );
        return atat::norm2( parent->origin->pos - (*i_bond)->origin->pos -
                            parent->structure->cell * (*i_translation) );
      }
      atat::rVector3d& vector( atat::rVector3d &_hold )
      {
        _hold = (*i_bond)->origin->pos - parent->origin->pos ;
        if ( *i_do_translate )
          _hold += parent->structure->cell * (*i_translation);
        return _hold;
      }
      types::t_real scalar_product( const const_iterator &_b ) const
      {
        atat::rVector3d a, b;
        if ( *i_do_translate )
          a =   (*i_bond)->origin->pos - parent->origin->pos 
              + parent->structure->cell * (*i_translation);
        else
          a = (*i_bond)->origin->pos - parent->origin->pos;
        if ( *_b.i_do_translate )
          b =   (*_b.i_bond)->origin->pos - _b.parent->origin->pos 
              + _b.parent->structure->cell * (*_b.i_translation);
        else
          b = (*_b.i_bond)->origin->pos - _b.parent->origin->pos;
        return a * b;
      }
      types::t_unsigned kind() const
        { return parent->bond_kind( *(*i_bond) ); }
      Ising_CE::Atom& get_origin()
        { return ((*i_bond)->get_origin()); }
      void translate( atat::rVector3d &_v, const atat::rMatrix3d &_cell )
        { if( *i_do_translate ) _v += _cell * ( *i_translation ); }
      void translate( atat::rVector3d &_v )
        { if( *i_do_translate ) _v += parent->structure->cell * ( *i_translation ); }
  }; // end of const_iterator definition

  // one for each site and type
  class Atomic_Functional 
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Vff::Atomic_Functional> ( Vff::Atomic_Functional& );
#endif
     const static types::t_real twos3;  // 2*sqrt(3)
     const static types::t_real one16;  // 1 / 16
     const static types::t_real s3o160; // sqrt(3) / 8
     const static types::t_real one640; // 1 / 640
     const static types::t_real three8; // 3 / 8
     const static types::t_real s33o8;  // 3 * sqrt(3) / 8
     const static types::t_real s33o16; // 3 * sqrt(3) / 16
     const static types::t_real thre16; // 3 / 8
     const static types::t_real thre32; // 3 / 32
     const static types::t_real s33128; // 3/128 *sqrt(3)
     const static types::t_real s33256; // 3/256 *sqrt(3)
     const static types::t_real no1280; 
     const static types::t_real no2560; 

    protected:
      std::string str;
      Ising_CE :: Structure *structure;
      types::t_unsigned site, type;
      std::vector< types::t_real > lengths; 
      std::vector< types::t_real > alphas; 
      std::vector< types::t_real > betas; 
      std::vector< types::t_real > gammas; 
      std::vector< types::t_real > sigmas; 
      
    public:
#ifdef _MPI
      Atomic_Functional   ( Ising_CE::Structure &_struct ) 
                        : structure(&_struct) {}  // should only be used by serialize
#endif
      Atomic_Functional   ( std::string _str, Ising_CE::Structure &_struct, 
                            types::t_unsigned _site, 
                            types::t_unsigned _type )
                        : str(_str), structure(&_struct), site(_site), type(_type) {}
      Atomic_Functional   ( const Atomic_Functional &_a )
                        : str(_a.str), structure(_a.structure), site(_a.site), type(_a.type),
                          lengths(_a.lengths), alphas(_a.alphas),
                          betas(_a.betas), gammas(_a.gammas), sigmas(_a.sigmas) {}
      
      void add_bond( const types::t_unsigned _typeB, const types::t_real _l,
                     const std::vector<types::t_real> &_i )
      {
        if ( lengths.size() < _typeB + 1)
          lengths.resize( _typeB+1, types::t_real(0) );
        if ( alphas.size() < (_typeB+1) * 5  )
          alphas.resize( (_typeB+1)*5, types::t_real(0) );
        lengths[_typeB] = _l;
        std::copy( _i.begin(), _i.end(), alphas.begin() + _typeB * 5 );
      }
      void add_angle( const types::t_unsigned _typeA,
                      const types::t_unsigned _typeC,
                      const types::t_real _gamma, const types::t_real _sigma, 
                      const std::vector<types::t_real> &_i )
      {
        types::t_unsigned offset = _typeA+_typeC;
        if ( gammas.size() < offset + 1  )
          gammas.resize( offset + 1, types::t_real(0) );
        if ( sigmas.size() < offset + 1  )
          sigmas.resize( offset + 1, types::t_real(0) );
        if ( betas.size() < (offset + 1) * 5 )
          betas.resize( (offset+1)*5, types::t_real(0) );
        gammas[offset] = _gamma;
        sigmas[offset] = _sigma;
        std::copy( _i.begin(), _i.end(), betas.begin() + offset * 5 );
      }
      
      types::t_real evaluate( const Atomic_Center &_center ) const;
      types::t_real evaluate_with_gradient( Atomic_Center &_center,
                                            const atat::rMatrix3d &_strain,
                                            atat::rMatrix3d &_stress,
                                            const atat::rMatrix3d &_K0 ) const;
      // computes the trace of the microscopic strain on an atomic center
      // structure0 and the atomic centers are expected to be related 
      types::t_real MicroStrain( const Atomic_Center &_center, 
                                 const Ising_CE::Structure &_str0 ) const;

      void print_out( std::ostream &stream ) const;
  }; 

  // actual vff functional
  class Functional : public function :: Base<>
  {
#ifdef _MPI
    friend bool mpi::BroadCast::serialize<Vff::Functional> ( Vff::Functional& );
#endif
    public:
      typedef types::t_real t_Type;
      typedef std::vector<t_Type> t_Container;
      typedef t_Container :: iterator iterator;
      typedef t_Container :: const_iterator const_iterator;

    protected:
      Ising_CE :: Structure &structure;
      Ising_CE :: Structure structure0; // needed for gradients
      types::t_real bond_cutoff;
      std::vector< Atomic_Center > centers;
      std::vector< Atomic_Functional > functionals;
      atat::rMatrix3d stress, strain;
      atat::rVector3d center_of_mass;
      
    public:
      Functional   ( Ising_CE :: Structure &_str )
                 : function::Base<>( 7 + _str.atoms.size() ),
                   structure(_str), structure0(_str), center_of_mass(0,0,0) {};
      Functional   ( const Vff::Functional &_c )
                 : function::Base<>( _c ), structure( _c.structure ),
                   structure0( _c.structure0 ), bond_cutoff( _c.bond_cutoff ),
                   centers( _c.centers ), functionals( _c.functionals ),
                   center_of_mass(_c.center_of_mass) {}
      ~Functional() {}

      bool Load( const TiXmlElement &_element );
      types::t_real evaluate(); // unpacks function::Base::variables, then calls energy
      types::t_real energy() const; 
      template< typename t_grad_iterator>
        void evaluate_gradient( t_grad_iterator const &_i_grad )
          { evaluate_with_gradient( _i_grad ); }
      void evaluate_gradient( t_Type * const _i_grad )
        { evaluate_with_gradient<t_Type*>( _i_grad ); }  
      template< typename t_grad_iterator>
        t_Type evaluate_with_gradient( t_grad_iterator const &_i_grad );
      t_Type evaluate_with_gradient( t_Type * const _i_grad )
        { return evaluate_with_gradient<t_Type*>( _i_grad ); }  
      t_Type evaluate_one_gradient( types::t_unsigned _pos) {return 0;}; // only meaningfull within VA
      bool init();
      void print_escan_input( const std::string &_f = "atom.config") const;
      bool construct_centers();
      bool initialize_centers();
      void print_out( std::ostream &stream ) const;
      const atat::rMatrix3d& get_stress() const
        { return stress; }

    protected:
      void unpack_variables(atat::rMatrix3d& strain);
      void pack_variables(const atat::rMatrix3d& _strain);
      template< typename t_grad_iterator>
      void pack_gradients(const atat::rMatrix3d& _stress, t_grad_iterator const &_grad) const;
  };

  template< typename t_grad_iterator>
  void Functional :: pack_gradients(const atat::rMatrix3d& _stress, 
                                    t_grad_iterator const &_grad) const
  {
    t_grad_iterator i_grad(_grad);

    // first, external stuff
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XX) )
      *i_grad = _stress(0,0), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_YY) ) 
      *i_grad = _stress(1,1), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_ZZ) ) 
      *i_grad = _stress(2,2), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XY) ) 
      *i_grad = 0.5 * (_stress(0,1) + _stress(1,0)), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_XZ) ) 
      *i_grad = 0.5 * (_stress(0,2) + _stress(2,0)), ++i_grad;
    if ( not (structure.freeze & Ising_CE::Structure::FREEZE_YZ) ) 
      *i_grad = 0.5 * (_stress(1,2) + _stress(2,1)), ++i_grad;

    // then atomic position stuff
    std::vector<Atomic_Center> :: const_iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: const_iterator i_end = centers.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom0 = structure0.atoms.begin();
    i_center = centers.begin();
    for (; i_center != i_end; ++i_center, ++i_atom0)
    {
      const atat::rVector3d& gradient = i_center->get_gradient();
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_X) ) 
        *i_grad = gradient[0], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Y) ) 
        *i_grad = gradient[1], ++i_grad;
      if ( not (i_atom0->freeze & Ising_CE::Atom::FREEZE_Z) ) 
        *i_grad = gradient[2], ++i_grad;
    }
  }

  template< typename t_grad_iterator>
  types::t_real Functional :: evaluate_with_gradient( t_grad_iterator const &_i_grad )
  {
    t_Type energy = 0;
    std::for_each( centers.begin(), centers.end(), std::mem_fun_ref(&Atomic_Center::reset_gradient) );

    // unpacks variables into vff atomic_center and strain format
    unpack_variables(strain);

    // computes K0
    atat::rMatrix3d K0 = (!(~strain));

    // computes energy and gradient
    std::vector<Atomic_Center> :: iterator i_center = centers.begin();
    std::vector<Atomic_Center> :: iterator i_end = centers.end();
    stress.zero();
    for (; i_center != i_end; ++i_center)
      energy += functionals[i_center->kind()].evaluate_with_gradient( *i_center, strain, stress, K0 );

    // now repacks into function::Base format
    pack_gradients(stress, _i_grad);

    return energy;
  }

} // namespace vff 

#endif // _VFF_FUNCTIONAL_H_
