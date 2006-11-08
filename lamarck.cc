#include "lamarck.h"

#include <lamarck/convex_hull.h>
#include <lamarck/one_point_hull.h>

#include <limits.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <math.h>

#include <eo/utils/eoRNG.h>

using opt::NO_MINIMIZER;
using opt::WANG_MINIMIZER;
using opt::PHYSICAL_MINIMIZER;
using opt::LINEAR_MINIMIZER;
using opt::SA_MINIMIZER;

namespace LaDa 
{
  Lamarck :: ~Lamarck()
  {
    if (convex_hull)
      delete convex_hull;
    convex_hull = NULL;
    if ( minimizers.size() )
    {
      std::vector< opt::Minimize_Base<t_GA_Functional> * > :: iterator i_minimizer = minimizers.begin();
      std::vector< opt::Minimize_Base<t_GA_Functional> * > :: iterator i_end = minimizers.end();
      for( ; i_minimizer != i_end; ++i_minimizer )
        delete *i_minimizer;
    }
  };

  bool Lamarck :: Load(const std::string &_filename) 
  {
    filename = _filename;
    TiXmlDocument doc( filename.c_str() );
    
    if  ( !doc.LoadFile() )
    {
      std::cerr << "error while opening input file " << filename << std::endl
                << doc.ErrorDesc() << std::endl; 
      return false;
    }

    TiXmlHandle docHandle( &doc );
    if ( not Load ( docHandle ) )
    {
      std::cerr << "error while  loading Lamarck parameters from " << filename << std::endl
                << " tinyxml: " << doc.ErrorDesc() << std::endl;
      return false;
    }

    // completes cluster list with equivalent clusters
    add_equivalent_clusters();
    
    read_CH();
    
    generate_functional(structure, &functional);
    
    // fitness pointers  -- note that functiona->variables will be set
    // to Individual::variables at each evaluation
    fitness.set_baseline( convex_hull );
    fitness.set_quantity( &functional ); 
    functional.destroy_variables(); 


    return true;
  }

  bool Lamarck :: read_CH()
  {
    { // first, creates ch
      if ( convex_hull )
        delete convex_hull;
      if ( is_one_point_hull )
        convex_hull = new VA_CE :: One_Point_Hull();
      else
      {
        convex_hull = new VA_CE :: Convex_Hull();
        init_convex_hull();
      }
    }
    { // then attempts to ch file
      TiXmlDocument doc( xml_filename.c_str() );
      
      if  ( !doc.LoadFile() ) // no convex hull to read...
        return false;
    
      TiXmlHandle handle( &doc );
      TiXmlElement *child = handle.FirstChild( "LaDa" ).FirstChild( "ConvexHull" ).Element();
      if ( not child )
        return false;
      
  
      if ( not convex_hull )
      {
        std::cerr << "Error while creating convex hull" << std::endl;
        return false;
      }
      
      return convex_hull->Load(child, *axes);
    } // end of attempt to read ch file
  }


  // loads input from filename file
  bool Lamarck :: Load( TiXmlHandle &handle )
  {
    TiXmlElement *child;
    rVector3d vec;

    // clusters, lattice, harmonics ....
    if ( not Functional_Builder :: Load (handle ) )
      return false;
     
    xml_filename = "convex_hull.xml";
    child = handle.FirstChild( "LaDa" ).FirstChild( "Filename" ).Element();
    if ( child and child->Attribute("xml") )
      xml_filename = child->Attribute("xml");
      
    // reads structure from input
    child = handle.FirstChild( "LaDa" ).FirstChild( "Structure" ).Element();
    if ( not child )
    {
      std::cerr << "Could not find Structure in " << filename << std::endl;
      return false;
    }
    structure.Load(child, *axes);

    child = handle.FirstChild( "LaDa" ).FirstChild( "GA" ).FirstChild("OnePointHull").Element();
    if ( child )
      is_one_point_hull = true;

    return true;
  }

  void Lamarck :: print_xml()
  {
    TiXmlDocument doc;
    TiXmlElement *LaDa_node = new TiXmlElement("LaDa");
    TiXmlElement *convex_hull_node = new TiXmlElement("ConvexHull");

    doc.SetTabSize(1);
    doc.LinkEndChild( new TiXmlDeclaration("1.0", "", "") );
    
    doc.LinkEndChild( LaDa_node );
    LaDa_node->LinkEndChild( convex_hull_node ); 
    convex_hull->print_xml( convex_hull_node, *axes );

    doc.SaveFile(xml_filename.c_str());
  }


  // adds endpoints to convex hull -- cluster list should be complete
  // at this point!!
  // single site CEs only
  void Lamarck :: init_convex_hull()
  {
    Ising_CE::Structure struc; 
    t_VA_Functional func;

    // creates structure...
    struc.cell = lattice->cell;
    struc.atoms.clear();
    struc.atoms.push_back( Ising_CE::Atom( lattice->atom_pos(0), -1.0 ) );

    // and functional - note that there is no CS 
    // (in fact CS may segfault here )
    generate_functional(struc, &func);
    

    // adds first endpoint
    {
      *(func.begin()) = -1.0;
      VA_CE::Breaking_Point bp(func.get_Obj1()->evaluate(), struc);
      static_cast<VA_CE::Convex_Hull*>(convex_hull)->force_add( bp );
    }

    // adds second endpoint
    {
      *(func.begin()) = 1.0;
      struc.atoms.clear();
      struc.atoms.push_back( Ising_CE::Atom( lattice->atom_pos(0), 1.0 ) );
      VA_CE::Breaking_Point bp(func.get_Obj1()->evaluate(), struc);
      static_cast<VA_CE::Convex_Hull*>(convex_hull)->force_add( bp );
    }

    // clean-up
    func.destroy_variables();
    delete func.get_Obj1();
    delete func.get_Obj2();
  }



  // adds a minimizer in list 
  t_unsigned Lamarck :: add_minimizer( t_unsigned _type, t_unsigned _n)
  {
    opt::Minimize_Base<t_GA_Functional> *minimizer;

    switch( _type )
    {
      case NO_MINIMIZER: 
        minimizer = new opt::Minimize_Base<t_GA_Functional>;
        break;
      case SA_MINIMIZER: 
        minimizer = new opt::Minimize_Linear<t_GA_Functional>;
        static_cast< opt::Minimize_Linear<t_GA_Functional>* >(minimizer)->simulated_annealing = true;
        static_cast< opt::Minimize_Linear<t_GA_Functional>* >(minimizer)->max_calls = _n;
        break;
      case LINEAR_MINIMIZER: 
        minimizer = new opt::Minimize_Linear<t_GA_Functional>;
        static_cast< opt::Minimize_Linear<t_GA_Functional>* >(minimizer)->simulated_annealing = false;
        static_cast< opt::Minimize_Linear<t_GA_Functional>* >(minimizer)->max_calls = _n;
        break;
      default:
        std::cerr << "Unknown minimizer type in LaDa :: MinimizerOp" << std::endl;
        exit(-1);
        break;
    }

    minimizers.push_back( minimizer );
    return minimizers.size() - 1;
  };

  bool Lamarck :: minimize( const t_Individual &_indiv, const t_unsigned &_minimizer ) 
  {
    if ( _minimizer > minimizers.size() ) 
      return false;

    fitness.set_variables( _indiv.get_variables() );
    return (minimizers[_minimizer])->minimize(fitness);
  }
  t_real Lamarck :: evaluate( t_Individual &_indiv ) 
  {
    ++EvalCounter;
    fitness.set_variables( _indiv.get_variables() );
    t_real result = ( fitness.get_quantity() )->evaluate();
    structure.set_atom_types( *_indiv.get_variables() );
    if ( convex_hull->add_structure(result, structure) )
      _indiv.invalidate_baseline();
    return result;
  }
  void Lamarck :: add_to_convex_hull( const t_Individual &_indiv ) 
  {
    structure.set_atom_types( *_indiv.get_variables() );
    if ( convex_hull->add_structure( _indiv.get_quantity(), structure) )
      _indiv.invalidate_baseline();
  }

  void Lamarck :: print_xmgrace( std::ofstream &_f, bool _print_ch )
  {
    std::string special_char = ""; 
    if ( not _print_ch ) 
      special_char = "? ";
    _f << " # " << special_char << " evaluation calls: " << EvalCounter << std::endl;
    _f << " # " << special_char << "polynomial calls: " 
       << VA_CE::Polynome::nb_eval << " "
       << VA_CE::Polynome::nb_eval_grad << " "
       << VA_CE::Polynome::nb_eval_with_grad << std::endl;
    _f << " # " << special_char << "bad guess " << opt::Minimize_Linear<t_GA_Functional> :: bad_guess 
       << "   good guess " << opt::Minimize_Linear<t_GA_Functional> :: good_guess 
       << std::endl
       << " # " << special_char << "poleval calls " << opt::Minimize_Linear<t_GA_Functional> :: nb_evals
       << "   polgrad calls " << opt::Minimize_Linear<t_GA_Functional> :: nb_grad_evals
       << std::endl;
 
    if( _print_ch )
      convex_hull->print_out(_f, VA_CE::Convex_Hull::PRINT_XMGRACE);
  }

  bool Lamarck :: Krossover( t_Individual  &_offspring, const t_Individual &_parent)
  {
    typedef std::complex<t_Individual :: t_Type> t_complex;
    typedef std::vector< t_complex > t_k_type;
    const t_complex imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    const std::vector<rVector3d> &k_vecs = get_kvectors(_offspring);
    typedef std::vector< t_complex > t_k_type;
    t_k_type k_offspring( k_vecs.size(), t_complex(0) );
    t_k_type k_parent( k_vecs.size(), t_complex(0) );
    
    // first, FTs parent and offspring
    std::vector<rVector3d> :: const_iterator i_kvec = k_vecs.begin();
    t_k_type :: iterator i_val = k_offspring.begin();
    t_k_type :: iterator i_val_end = k_offspring.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom_begin = structure.atoms.begin();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom_end = structure.atoms.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom;
    t_Individual :: const_iterator i_spin_begin = _offspring.begin();
    t_Individual :: const_iterator i_spin;
    std::cout << std::setprecision(2);
    for (t_int i=0; i < 2; ++i)
    {
      for ( ; i_val != i_val_end; ++i_val, ++i_kvec)
      {
        i_atom = i_atom_begin;
        i_spin = i_spin_begin;
        for(; i_atom != i_atom_end; ++i_atom, ++i_spin )
        {
          *i_val +=    exp( imath * ( i_atom->pos[0] * (*i_kvec)[0] +
                                      i_atom->pos[1] * (*i_kvec)[1] +
                                      i_atom->pos[2] * (*i_kvec)[2] ) )
                     * (*i_spin);
        }
      }
      i_val = k_parent.begin();    // FT _parent next
      i_val_end = k_parent.end();
      i_spin_begin = _parent.begin();
      i_kvec = k_vecs.begin();
    }

    // then does crossover
    t_k_type :: const_iterator i_cnst = k_parent.begin();
    i_val = k_offspring.begin();
    i_val_end = k_offspring.end();
    for ( ; i_val != i_val_end; ++i_val, ++i_cnst)
      if ( rng.flip() )
        *i_val = *i_cnst;

    // Then FT back to r space, while making sure values are +/-1
    std::vector<rVector3d> :: const_iterator i_kvec_begin = k_vecs.begin();
    t_k_type :: iterator i_val_begin = k_offspring.begin();
    t_Individual :: iterator i_var = _offspring.begin();
    i_atom = i_atom_begin;
    for ( ; i_atom != i_atom_end; ++i_var, ++i_atom)
    {
      t_complex store(0);
      i_kvec = i_kvec_begin;
      i_val = i_val_begin;
      for(; i_val != i_val_end; ++i_kvec, ++i_val )
      {
        store +=    exp( -imath * ( i_atom->pos[0] * (*i_kvec)[0] +
                                    i_atom->pos[1] * (*i_kvec)[1] +
                                    i_atom->pos[2] * (*i_kvec)[2] ) )
                  * (*i_val);
      }
      *i_var = ( std::real( store ) > 0 ) ? 1.0 : -1.0;
    }

    return true; // offspring has changed!
  }

  void Lamarck :: fourrier_transform( const t_Individual &_indiv, 
                                      std::vector< std::complex<t_Individual :: t_Type> > &_fourrier )
  {
    typedef std::complex<t_Individual :: t_Type> t_complex;
    typedef std::vector< t_complex > t_k_type;
    const t_complex imath(0, -2*3.1415926535897932384626433832795028841971693993751058208);
    const std::vector<rVector3d> &k_vecs = get_kvectors(_indiv);
    std::vector<rVector3d> :: const_iterator i_kvec = k_vecs.begin();
    t_k_type :: iterator i_val = _fourrier.begin();
    t_k_type :: iterator i_val_end = _fourrier.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom_begin = structure.atoms.begin();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom_end = structure.atoms.end();
    std::vector<Ising_CE::Atom> :: const_iterator i_atom;
    t_Individual :: const_iterator i_spin_begin = _indiv.begin();
    t_Individual :: const_iterator i_spin;

    _fourrier.resize( k_vecs.size() );

    for ( ; i_val != i_val_end; ++i_val, ++i_kvec)
    {
      i_atom = i_atom_begin;
      i_spin = i_spin_begin;
      *i_val = t_complex(0);
      for(; i_atom != i_atom_end; ++i_atom, ++i_spin )
      {
        *i_val +=    exp( imath * ( i_atom->pos[0] * (*i_kvec)[0] +
                                    i_atom->pos[1] * (*i_kvec)[1] +
                                    i_atom->pos[2] * (*i_kvec)[2] ) )
                   * (*i_spin);
      }
    }
  }
  
} // namespace LaDa
