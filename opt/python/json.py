""" Module to decorate properties with json transcripters. """

def unit(unit):
  """ Creates JSON transfer functions wich remove/add units. """
  def to_json(object):
    """ Removes unit from object. """
    return float(object.magnitude)
  def from_json(object):
    """ Adds unit to object. """
    return object * unit
  def add_json_transfer_functions(function):
     """ Adds json transfer functions to an object. """
     function.to_json = to_json
     function.from_json = from_json
     return function
  return add_json_transfer_functions

def array(type):
  """ Creates JSON transfer functions wich transforms numpy arrays to list. """
  def to_json(object):
    """ Transforms array to list. """
    return object.tolist()
  def from_json(object):
    """ Transforms list to array. """
    from numpy import array
    return array(object, dtype=type)
  def add_json_transfer_functions(function):
     """ Adds json transfer functions to an object. """
     function.to_json = to_json
     function.from_json = from_json
     return function
  return add_json_transfer_functions

def array_with_unit(type, unit):
  """ Creates JSON transfer functions wich transforms numpy arrays to list. """
  def to_json(object):
    """ Transforms array to list. """
    return object.magnitude.tolist()
  def from_json(object):
    """ Transforms list to array. """
    from numpy import array
    return array(object, dtype=type) * unit
  def add_json_transfer_functions(function):
     """ Adds json transfer functions to an object. """
     function.to_json = to_json
     function.from_json = from_json
     return function
  return add_json_transfer_functions


def pickled(function):
  """ Adds JSON transfer functions which work through a pickle. """
  from pickle import dumps, loads
  def to_json(object):
    """ Transform to pickle. """
    return dumps(object)
  def from_json(object):
    """ Transform to pickle. """
    return loads(str(object))
  function.to_json = to_json
  function.from_json = from_json
  return function
