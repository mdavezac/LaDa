""" Class to handle creation/destruction of temporary attributes. 
   
    If you want to create a temporary attribute in self within a given scope:
    >>> with Tempdir(self, "attrname", attrvalue) as self:
    >>>   ...
"""
class Tempattr:
  """ Works with "with" statement to create/destroy a temporary attribute.  """
  def __init__(self, instance, name, value):
    self.instance = instance
    self.name = name
    self.value = value


  def __enter__(self):
    """ Creates temporary attribute """

    self.hasvalue = hasattr(self.instance, self.name)
    if self.hasvalue: self.oldvalue = getattr(self.instance, self.name)
    setattr(self.instance, self.name, self.value)
    return self.instance

  def __exit__(self, type, value, traceback):
    """ Deletes temporary attribute """
    if self.hasvalue: setattr(self.instance, self.name, self.oldvalue)
    else: delattr(self.instance, self.name)
