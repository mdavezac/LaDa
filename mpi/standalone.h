//
//  Version: $Id$
//
#ifndef _MPI_STANDALONE_H_
#define _MPI_STANDALONE_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <opt/types.h>

#include "comm.h"

namespace mpi
{

  //! Sends a template object from current process to \a _bull.
  template< class T_OBJECT >
  void send_template_object( types::t_unsigned _bull,
                             const T_OBJECT &_o, MPI::Intracomm *_comm);
  //! Receives a template object from \a _bull.
  template< class T_OBJECT >
  void receive_template_object( types::t_unsigned _bull,
                                T_OBJECT &_o, MPI::Intracomm *_comm);
  //! Sends a (non-template) object to bull \a _bull.
  template< class T_OBJECT >
  void send_object( types::t_unsigned _bull,
                    const T_OBJECT _o, MPI::Intracomm *_comm);
  //! Receives a (non-template) object from bull \a _bull.
  template< class T_OBJECT >
  void receive_object( types::t_unsigned _bull,
                       const T_OBJECT &_o, MPI::Intracomm *_comm);
  //! Sends a range (must be iterators or pointers).
  template< class T_ITERATOR >
  void send_range( types::t_unsigned _bull, const T_ITERATOR _first,
                   const T_ITERATOR  _last, MPI::Intracomm *_comm);
  //! Receive a range (must be iterators or pointers).
  template< class T_ITERATOR >
  void receive_range( types::t_unsigned _bull, T_ITERATOR _first,
                      T_ITERATOR  _last, MPI::Intracomm *_comm);
  //! Broadcasts a template object from \a _root.
  template< class T_OBJECT >  
  void bcast_template_object( types::t_unsigned _root,
                              T_OBJECT &_object, MPI::Intracomm *_comm );
  //! Broadcasts a template constant object from \a _root.
  template< class T_OBJECT >  
  void const_bcast_template_object( types::t_unsigned _root,
                                    const T_OBJECT &_object, MPI::Intracomm *_comm );
  template< class T_ITERATOR >  
    void bcast_range( types::t_unsigned _root, T_ITERATOR _first,
                      T_ITERATOR _last, MPI::Intracomm *_comm );
  //! Broadcasts a range of constant iterators/pointers.
  template< class T_ITERATOR >  
    void const_bcast_range( types::t_unsigned _root, const T_ITERATOR _first,
                            const T_ITERATOR _last, MPI::Intracomm *_comm );
  //! Broadcasts a template object from \a _root.
  template< class T_OBJECT >  
  void bcast_object( types::t_unsigned _root, T_OBJECT &_object, MPI::Intracomm *_comm );
  //! Broadcasts a constant template object from \a _root.
  template< class T_OBJECT >  
  void const_bcast_object( types::t_unsigned _root, const T_OBJECT &_object, 
                           MPI::Intracomm *_comm );

    template< class T_ITERATOR >
      void send_range( types::t_unsigned _bull, const T_ITERATOR _first,
                       const T_ITERATOR _last, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc.serialize( _first, _last );
        bc.allocate_buffers();
        bc.serialize( _first, _last );
        bc.send_ptp( _bull );
      }

    template< class T_ITERATOR >
      void receive_range( types::t_unsigned _bull, T_ITERATOR _first,
                          T_ITERATOR _last, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc.serialize( _first, _last );
        bc.allocate_buffers();
        bc.serialize( _first, _last );
        bc.receive_ptp( _bull );
        bc.serialize( _first, _last );
      }

    template< class T_OBJECT >  
      void bcast_template_object( types::t_unsigned _root, T_OBJECT &_object,
                                  MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        _object.serialize( bc );
        bc.allocate_buffers();
        _object.serialize( bc ); 
        bc( _root );
        _object.serialize( bc );
      }
    template< class T_OBJECT >  
      void const_bcast_template_object( types::t_unsigned _root, const T_OBJECT &_object,
                                        MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        _object.serialize( bc );
        bc.allocate_buffers();
        _object.serialize( bc ); 
        bc( _root );
      }
    template< class T_ITERATOR >  
      void bcast_range( types::t_unsigned _root, T_ITERATOR _first,
                        T_ITERATOR _last, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc.serialize( _first, _last );
        bc.allocate_buffers();
        bc.serialize( _first, _last );
        bc( _root );
        bc.serialize( _first, _last );
      }
    template< class T_ITERATOR >  
      void const_bcast_range( types::t_unsigned _root, const T_ITERATOR _first,
                              const T_ITERATOR _last, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc.serialize( _first, _last );
        bc.allocate_buffers();
        bc.serialize( _first, _last );
        bc( _root );
      }
    template< class T_OBJECT >  
      void bcast_object( types::t_unsigned _root, T_OBJECT &_object, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc << _object << BroadCast::allocate
           << _object; 
        bc( _root ); 
        bc << _object;
      }
    template< class T_OBJECT >  
      void const_bcast_object( types::t_unsigned _root, const T_OBJECT &_object, 
                               MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc << _object << BroadCast::allocate
           << _object;
        bc( _root ); 
      }



    template< class T_OBJECT >  
      void receive_template_object( types::t_unsigned _bull, T_OBJECT &_object,
                                    MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        _object.serialize( bc );
        bc.allocate_buffers();
        _object.serialize( bc );
        bc.receive_ptp( _bull );
        _object.serialize( bc );
      }

    template< class T_OBJECT >  
      void send_template_object( types::t_unsigned _bull, const T_OBJECT &_object,
                                 MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        std::cout << "sending template object " << _object << std::endl;
        BroadCast bc( _comm );
        _object.serialize( bc );
        bc.allocate_buffers();
        _object.serialize( bc );
        bc.send_ptp( _bull );
      }
    template< class T_OBJECT >  
      void receive_object( types::t_unsigned _bull,
                           T_OBJECT &_object, MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc << _object << BroadCast::allocate << _object;
        bc.receive_ptp( _bull );
        bc << _object;
      }

    template< class T_OBJECT >  
      void send_object( types::t_unsigned _bull, const T_OBJECT _object,
                        MPI::Intracomm *_comm )
      {
        if( _comm->Get_size() < 2 ) return;
        BroadCast bc( _comm );
        bc << _object << BroadCast::allocate << _object;
        bc.send_ptp( _bull );
      }

}

#endif
