//
//  Version: $Id$
//
namespace LaDa
{
  namespace GA
  {
    namespace PureLayers
    {
  
      template< class T_OBJECT >
        void Concentration :: operator()( T_OBJECT &_obj )
        {
          // computes concentrations first
          get( _obj );
          if( not single_c )
          {
            x0 = x; y0 = y;
            if ( rng.flip() or not can_inverse(x) ) x0 = get_x(y0);
            else                                    y0 = get_y(x0);
          }
    
          // compute number of changes to make per site type
          types::t_real to_change[2];
          to_change[0] = (types::t_real) N * ( x0 - x );
          to_change[1] = (types::t_real) N * ( y0 - y );
          if(     math::gt(to_change[0],  -1.0)  and math::le(to_change[0], 1.0 )
              and math::gt(to_change[1],  -1.0)  and math::le(to_change[1], 1.0 ) ) return;
    
          __ASSERT( sites.size() != _obj.Container().size(), "Inequivalent sizes.\n" )
    
          // Puts the positions which can be changed into a list
          std::vector<types::t_unsigned> pos[2];
          typedef typename T_OBJECT :: t_Container :: const_iterator const_iterator;
          std::vector<bool> :: const_iterator i_site = sites.begin();
          const_iterator i_bit = _obj.Container().begin();
          const_iterator i_bit_end = _obj.Container().end();
          types::t_unsigned a = 0, b = 0;
          for(types::t_unsigned i=0; i_bit != i_bit_end; ++i_bit, ++i, ++i_site)
          {
            if( *i_site ) ++a;
            else ++b;
            types::t_unsigned site_index = *i_site ? 0: 1;
            if( to_change[site_index] > 0.0 and BitString::spin_down(*i_bit)  )
              pos[site_index].push_back( i );
            else if( to_change[site_index] < 0.0 and BitString::spin_up(*i_bit)  )
              pos[site_index].push_back( i );
          }
    
          // shuffle position lists
          types::t_unsigned (*ptr_to_rng)(const types::t_unsigned& )
              = &eo::random<types::t_unsigned>;
          for( types::t_unsigned i=0; i < 2; ++i )
          {
            if( math::gt(to_change[i],  -1.0) and math::le(to_change[i], 1.0 ) ) continue;
            std::vector<types::t_unsigned> :: iterator i_pos = pos[i].begin();
            std::vector<types::t_unsigned> :: iterator i_pos_end = pos[i].end();
            std::random_shuffle(i_pos, i_pos_end, ptr_to_rng );
            i_pos = pos[i].begin();
            for(; i_pos != i_pos_end; ++i_pos)
            {
              _obj.bitstring[*i_pos] = _obj.bitstring[*i_pos] > 0e0 ? -1e0: 1e0;
              ( to_change[i] > 0 ) ? to_change[i] -= 2: to_change[i] += 2;
            
              // Fuzzy math at this point could create infinit loop
              if( math::gt(to_change[i],  -1.0) and math::leq(to_change[i], 1.0 ) ) break;
            }
            __DOASSERT( math::leq(to_change[i],  -1.0) or math::gt(to_change[i], 1.0 ),
                           "Concentration could not be set\n"
                        << "Incompatibility between required x/y and frozen atoms?\n"; )
          }
        }

      template< class T_OBJECT >
        void Concentration :: get( const T_OBJECT &_obj )
        {
          __ASSERT( sites.size() != _obj.bitstring.size(),
                    "Unequal sizes.\n" )
    
          // computes concentrations first
          typedef typename T_OBJECT :: t_Container :: const_iterator t_cit;
          t_cit i_bit = _obj.bitstring.begin();
          t_cit i_bit_end = _obj.bitstring.end();
          std::vector<bool>::const_iterator i_site = sites.begin();
          types::t_int concx = 0, concy = 0;
          for(; i_bit != i_bit_end; ++i_bit, ++i_site )
            if( *i_site ) *i_bit > 0 ? ++concx: --concx;
            else          *i_bit > 0 ? ++concy: --concy;
    
          // add frozen bits
          concx += Nfreeze_x;
          concy += Nfreeze_y;
    
          // finally normalizes
          x = (types::t_real) concx / (types::t_real) N;
          y = (types::t_real) concy / (types::t_real) N;
        }
  
      } // namespace PureLayers
  
  
    } // namespace GA


} // namespace LaDa
