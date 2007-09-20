//
//  Version: $Id$
//
#ifndef _RANKING_H_
#define _RANKING_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <tinyxml/tinyxml.h>

#include <opt/types.h>
#ifdef _MPI
#include <mpi/mpi_object.h>
#endif

namespace Ranking
{
  template<class T_GATRAITS> 
    class Base : public eoUF< typename T_GATRAITS::t_Population&, void> {};

  template<class T_SHARING>
  class Niching : Base
  {
    public:
      typedef T_SHARING t_Sharing;
      typedef typename t_Sharing :: t_GAOpTraits t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef typename t_IndivTraits :: t_Fitness t_Fitness;
      typedef typename t_Fitness :: t_Quantity t_FitnessQuantity;
      typedef typename t_IndivTraits :: t_Population t_Population;
      typedef typename t_IndivTraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;
      typedef std::vector< t_ScalarFitnessQuantity >  t_Column;

    protected:
      t_Sharing *sharing;
      t_Column sums;

    public:
  
      explicit Niching( t_Sharing *_share ) : t_Base(), sharing(_share) {}
      ~Niching() { if (sharing) delete sharing; sharing = NULL; } 
  
      void operator()(t_Population& _pop)
      {
        map( _pop );

        typename t_Population :: iterator i_indiv = _pop.begin();
        typename t_Population :: iterator i_end = _pop.end();
        typename t_Column :: const_iterator i_sum = sums.begin();
        for(; i_indiv != i_end; ++i_indiv );
          i_indiv->set_fitness( (t_ScalarFitnessQuantity&) i_indiv.fitness() / *i_sum );
      }

    protected:
      void map(const t_Population )
      {
        sums.clear(); sums.resize( _pop.size(), t_ScalarFitnessQuantity(0) );
        typename t_Population :: iterator i_indiv = _pop.begin();
        typename t_Population :: iterator i_end = _pop.end();
        for(unsigned i=0; i_indiv != i_end; ++i_indiv, ++i );
        {
          i_2indiv = i_indiv;  
          for(unsigned j=i; i_2indiv != i_end; ++i_2indiv, ++j );
          {
            sums[i] += (*share)(*i_indiv, *i_2indiv);
            sums[j] += (*share)(*i_indiv, *i_2indiv);
          }
        }
      }
  
    private :
  
    eoDominanceMap<EOT>& dominanceMap;
  };

}

namespace Distance
{
  template<class T_GAOPTRAITS>
  class Hamming
  {
    public:
      typedef T_GAOPTRAITS t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits::t_Object t_Object;
      typedef typename t_IndivTraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

      t_FitnessQuantity operator()( const t_Individual &_i1, const t_Individual &_i2) const
      {
        if ( _i1.Object().Container().size() != _i1.Object().Container().size() )
          throw std::runtime_error("individuals of different size in Distance::Hamming!\n");
        typename t_Object::const_iterator i_bit2 = _i2.Object().begin();
        typename t_Object::const_iterator i_bit1 = _i1.Object().begin();
        typename t_Object::const_iterator i_bit1_end = _i1.Object().end();
        t_ScalarFitnessQuantity result;
        for(; i_bit1 != i_bit1_end; ++i_bit1, ++i_bit2 )
          result += std::abs(   t_ScalarFitnessQuantity(*i_bit1) 
                              - t_ScalarFitnessQuantity(*i_bit2) );
        return result;
      }
  };
}

namespace Sharing
{
  template <class T_DISTANCE>
  class Triangular 
  {
    public:
      typedef T_DISTANCE t_Distance;
      typedef typename t_Distance :: t_GAOpTraits t_GAOpTraits;
    protected:
      typedef typename t_GAOpTraits :: t_Individual t_Individual;
      typedef typename t_GAOpTraits :: t_IndivTraits t_IndivTraits;
      typedef typename t_GAOpTraits :: t_Concentration t_Concentration;
      typedef typename t_IndivTraits :: t_Fitness t_Fitness;
      typedef typename t_IndivTraits :: t_ScalarFitness t_ScalarFitness;
      typedef typename t_ScalarFitness :: t_Quantity t_ScalarFitnessQuantity;

    protected:
      t_Distance *distance;
      types::t_unsigned alpha;

    public:
      explicit
        Triangular( t_Distance *_d ) : distance(_d), alpha(1)  {}
      ~Triangular() { if( distance ) delete distance; distance = NULL; }

      bool Load( const TiXmlElement &_node )
      {
        // First finds correct node
        std::string name = _node.Value();
        const TiXmlElement *parent;
        if( name == "Sharing" )
          parent = &_node;
        else
          parent = _node.FirstChildElement("Sharing");
        
        for(; parent; parent = parent->NextSiblingElement("Sharing") )
        {
          if ( not parent->Attribute("type") ) continue;
          name = parent->Attribute("type");
          if(    name == "Triangular"
              or name == "triangular" ) break;
        }

        if ( not parent ) return false;
        if ( parent->Attribute("alpha") )
        {
          int d = 0;
          parent->Attribute("alpha", d );
          if ( d < 0 ) 
          {
            std::cerr << "Invalid alpha = " << d 
                      << " in Sharing::Triangular " << std::endl;
            return false;
          } 
            
          alpha = types::t_unsigned( std::abs(d) );
        }
        if ( parent->Attribute("sigma") )
        {
          double d = 0;
          parent->Attribute("sigma", d );
          if ( d < 0 ) 
          {
            std::cerr << "Invalid sigma = " << d 
                      << " in Sharing::Triangular " << std::endl;
            return false;
          } 
            
          sigma = t_ScalarFitnessQuantity( d );
        }
      
        std::string print_out() const
        { 
          if ( not distance ) return "Unitialized Triangular Sharing"; 
          std::ostringstream sstr; 
          sstr << "Triangular Sharing, alpha=" << alpha << ", sigma=" << sigma
               << ", distance: " << distance->print_out();
          return sstr.str();
        }

        bool is_valid() const { return distance != NULL; }

        t_ScalarFitnessQuantity operator()(const t_Individual& _i1, const t_Individual& _i2 ) const
        {
          t_ScalarFitnessQuantity d = (*distance)( _i1, _i2 );
          if ( d > sigma ) return t_ScalarFitnessQuantity(0);
          if( alpha == -1 )
            return t_ScalarFitnessQuantity(1) - d / sigma;
          return std::pow( t_ScalarFitnessQuantity(1) - d / sigma, alpha );
        }
                                
      };
  }  
}



#endif // _RANKING_H_
