""" Class which keeps track of known candidates. """
__docformat__ = "restructuredtext en"
class History(object):
   """ Class which keeps track of known candidates. 
   
       The history is kept via a history file. 
   """
   def __init__(self, filepath):
     """ Initializes a history object. """
     self.filepath = filepath
     """ file where history is saved. """

   def remove_stale(self, comm=None):
     """ Remove possible stale lock. """
     from ..opt import LockFile
     LockFile(self.filepath).remove_stale(comm)

   def __call__(self, indiv, add=True):
     """ Returns true if individual is in history. """
     from pickle import load, dump
     from ..opt import acquire_lock
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       result = indiv in history
       if add and not result:
         history.append(history)
         with open(self.filepath, "r") as file: dump(history, file)
     return result

   def get(self, indiv):
     """ Returns individual if it exists in history. """
     from pickle import load
     from ..opt import acquire_lock
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       if indiv in history: return history[history.index(indiv)]

   def pop(self, indiv):
     """ Removes individual if it exists in history. """
     from pickle import load, dump
     from ..opt import acquire_lock
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       if indiv in history: 
         result = history.pop(history.index(indiv))
         with open(self.filepath, "r") as file: dump(history, file)
         return result

   def __repr__(self):
     """ Python representation of this object. """
     return "{0.__class__.__name__}({1})".format(self, repr(self.filepath))
