from ._crystal import AtomSet, AtomVec, AtomStr

def cast(self, kind, force=False):
  """ Casts an atom from one kind to another. 
  
      :param kind: The new kind to which to cast self. Can be "scalar", "list",
        or "set". If this is the same kind as self, then self is returned.
      :param force: Will not throw when downcasting to scalar kind encurs a
        loss of information.
  """
  from copy import deepcopy
  from .. import error

  assert kind in ["scalar", "list", "set"],\
         error.ValueError("Uknown kind {0}.".format(kind))
  if    (kind == "scalar" and self.__class__ == AtomStr)\
     or (kind == "list" and self.__class__ == AtomVec)\
     or (kind == "set" and self.__class__ == AtomSet): return self\

  if kind == "scalar":
    if not (force or len(self.type) < 2): 
      error.ValueError("Cannot downcast without loosing information.\n"\
                       "Use keyword argument ``force=True`` if you really "\
                       "really want to do this.")
    type = "" if len(self.type) == 0 else self.type.__iter__().next()
    result = AtomStr(*self.pos, type=type, site=self.site, freeze=self.freeze)
  elif kind == "list":
    type = [self.type] if self.kind == "scalar" else self.type
    result = AtomVec(*self.pos, type=type, site=self.site, freeze=self.freeze)
  elif kind == "set":
    type = set([self.type] if self.kind == "scalar" else self.type)
    result = AtomSet(*self.pos, type=type, site=self.site, freeze=self.freeze)

  result.__dict__.update( deepcopy(self.__dict__) )
  return result

AtomSet.cast = cast
AtomVec.cast = cast
AtomStr.cast = cast

def _repr(self):
  """ Represents an atom. """
  result = "Atom({0.pos[0]}, {0.pos[1]}, {0.pos[2]}".format(self)
  if self.kind == "scalar": result += ", {0}".format(repr(self.type))
  elif self.kind == "list": result += ", {0}".format(repr(list(self.type)))
  elif self.kind == "set": result += ", {0}, kind=\"set\"".format(repr(list(self.type)))
  if self.site != -1: result += ", site={0}".format(self.site)
  if self.freeze != 0: result += ", freeze={0}".format(self.freeze)
  result += ")"
  return result
AtomStr.__repr__ = _repr
AtomVec.__repr__ = _repr
AtomSet.__repr__ = _repr
