#ifndef LADA_PYTHON_RANDOM_ACCESS_LIST_ITERATOR_H
#define LADA_PYTHON_RANDOM_ACCESS_LIST_ITERATOR_H

#include <LaDaConfig.h>

#include <vector>

#include "object.h"

namespace LaDa 
{
  namespace python
  {
    //! \brief Thin wrapper around a python list implementing random access stl iterators.
    //! \details No reference count management. The point is that the list
    //!          should be owned and stable throughout the life of the iterator. 
    class RAList_iterator
    {
      public:
        //! \typedef type representing distance between iterators.
        typedef int difference_type;
        //! \typedef type representing returned by dereferencing.
        typedef PyObject* value_type;
        //! \typedef type representing pointer to referee.
        typedef std::vector<value_type>::iterator::pointer pointer;
        //! \typedef type representing reference to referee.
        typedef std::vector<value_type>::iterator::reference reference;
        //! \typedef type representing the random access iterator category.
        typedef std::vector<value_type>::iterator::iterator_category iterator_category;
        //! Initializes a random access iterator.
        RAList_iterator(PyObject* _in) : object_(_in), index_(PyList_GET_SIZE(_in)) {}
        //! Initializes a random access iterator.
        RAList_iterator(PyObject* _in, difference_type _index) : object_(_in), index_(_index) {}
        //! Copies a python reference.
        RAList_iterator(RAList_iterator const &_in) : object_(_in.object_), index_(_in.index_) {}
        //! Deref operator.
        const reference operator*() const
        {
          if(index_ >= (int)PyList_GET_SIZE(object_) or index_ < 0)
            LADA_PYTHROW(IndexError, "List iterator is out of range.");
          return PyList_GET_ITEM(object_, index_);
        }
        //! Deref operator.
        reference operator*()
        {
          if(index_ >= (int)PyList_GET_SIZE(object_) or index_ < 0)
            LADA_PYTHROW(IndexError, "List iterator is out of range.");
          return PyList_GET_ITEM(object_, index_);
        }
        //! Deref operator.
        const pointer operator->() const
        {
          if(index_ >= (int)PyList_GET_SIZE(object_) or index_ < 0)
            LADA_PYTHROW(IndexError, "List iterator is out of range.");
          return &PyList_GET_ITEM(object_, index_);
        }
        //! Deref operator.
        pointer operator->()
        {
          if(index_ >= (int)PyList_GET_SIZE(object_) or index_ < 0)
            LADA_PYTHROW(IndexError, "List iterator is out of range.");
          return &PyList_GET_ITEM(object_, index_);
        }
        //! Equality operator.
        bool operator==(RAList_iterator const &_in) const
          { return object_ == _in.object_ and index_ == _in.index_; }
        //! Inequality operator.
        bool operator!=(RAList_iterator const &_in) const
          { return object_ != _in.object_ or index_ != _in.index_; }
        //! Less than operator.
        bool operator<(RAList_iterator const &_in) const
          { return index_ < _in.index_; }
        //! Less than or equal operator.
        bool operator<=(RAList_iterator const &_in) const
          { return index_ <= _in.index_; }
        //! Greater than operator.
        bool operator>(RAList_iterator const &_in) const
          { return index_ > _in.index_; }
        //! Greater than or equal operator.
        bool operator>=(RAList_iterator const &_in) const
          { return index_ >= _in.index_; }
        //! Prefix Increment.
        //! Prefix Increment.
        //! Prefix Increment.
        RAList_iterator& operator++() { ++index_; return *this;}
        //! Postfix Increment.
        RAList_iterator operator++(int) { RAList_iterator tmp(*this); ++index_; return tmp;}
        //! Prefix Decrement.
        RAList_iterator& operator--() { --index_; return *this;}
        //! Postfix Decrement.
        RAList_iterator operator--(int) { RAList_iterator tmp(*this); --index_; return tmp;}
        //! Increments.
        void operator+=(difference_type _a) { index_ += _a; }
        //! Decrements.
        void operator-=(difference_type _a) { index_ += _a; }

        //! Increments.
        RAList_iterator operator+(difference_type _a) { return RAList_iterator(object_, index_ + _a); }
        //! Decrements.
        RAList_iterator operator-(difference_type _a) { return RAList_iterator(object_, index_ - _a); }
        //! Difference between pointers.
        difference_type operator-(RAList_iterator const &_a) { return index_ - _a.index_; }

      protected:
        //! Python to list.
        PyObject* object_;
        //! Current Item.
        difference_type index_;
    };
  }
}
#endif