""" Class which keeps track of known candidates. """
__docformat__ = "restructuredtext en"
class History(object):
   """ Class which keeps track of known candidates. 
   
       The history is kept via a history file. 
   """
   def __init__(self, filename):
     """ Initializes a history object. """
     self.filename = filename
     """ file where history is saved. """
 
   def remove_stale(self, comm=None):
     """ Remove possible stale lock. """
     from ..opt import LockFile
     LockFile(self.filename).remove_stale(comm)

   def __call__(self, indiv, add=True):
     """ Returns true if individual is in history. """
     from pickle import load, dump
     from ..opt import acquire_lock
     with acquire_lock(self.filename) as lock:
       with open(self.filename, "r") as file: history = load(file)
       result = indiv in history
       if add and not result:
         history.append(history)
         with open(self.filename, "r") as file: dump(history, file)
     return result

   def get(self, indiv):
     """ Returns individual if it exists in history. """
     from pickle import load
     from ..opt import acquire_lock
     with acquire_lock(self.filename) as lock:
       with open(self.filename, "r") as file: history = load(file)
       if indiv in history: return history[history.index(indiv)]

   def pop(self, indiv):
     """ Removes individual if it exists in history. """
     from pickle import load, dump
     from ..opt import acquire_lock
     with acquire_lock(self.filename) as lock:
       with open(self.filename, "r") as file: history = load(file)
       if indiv in history: 
         result = history.pop(history.index(indiv))
         with open(self.filename, "r") as file: dump(history, file)
         return result


