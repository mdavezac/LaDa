""" Class which keeps track of known candidates. """
__docformat__ = "restructuredtext en"

class History(object):
   """ Class which keeps track of known candidates. 
   
       The history is kept via a history file. 
   """
   HISTORYCAR = "HISTORYCAR"
   """ Default history filename. """

   def __init__(self, filename = None, directory=None):
     """ Initializes a history object. """
     from ..opt import RelativeDirectory

     self.filename = filename if filename != None else self.HISTORYCAR
     """ file where history is saved. """
     self._directory = RelativeDirectory(directory)
     """ Directory where history file is saved. """

   @property
   def directory(self):
     """ Directory where history file is saved. """
     return self._directory.path
   @directory.setter 
   def directory(self, value): self._directory.path = value

   @property
   def filepath(self):
     """ Path to history file. """
     from os.path import join
     return join(self.directory, self.filename)

   def remove_stale(self, comm=None):
     """ Remove possible stale lock. """
     from ..opt import LockFile
     LockFile(self.filepath).remove_stale(comm)

   def __call__(self, indiv, add=True):
     """ Returns true if individual is in history. """
     from pickle import load, dump
     from os.path import exists
     from ..opt import acquire_lock
     with acquire_lock(self.filepath) as lock:
       if not exists(self.filepath): history, result = [indiv], False
       else:
         with open(self.filepath, "r") as file: history = load(file)
         try: index = history.index(indiv)
         except ValueError:
           if not add: return False
           history.append(indiv)
           result = False
         else: # updates individual when possible.
           history[index].__dict__.update(indiv.__dict__)
           result = True
       with open(self.filepath, "w") as file: dump(history, file)
     return result

   def get(self, indiv):
     """ Returns individual if it exists in history. """
     from pickle import load
     from os.path import exists
     from ..opt import acquire_lock
     if not exists(self.filepath): return 
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       if indiv in history: return history[history.index(indiv)]

   def replace(self, indiv):
     """ Replaces an already existing indiviual. 


         Throws ValueError if the individual does not exist. 
     """ 
     from pickle import load
     from os.path import exists
     from ..opt import acquire_lock
     if not exists(self.filepath): return 
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       try: history[history.index(indiv)] = indiv
       except ValueError: raise ValueError("Individual not in history.")
       with open(self.filepath, "w") as file: dump(history, file)
       

   def pop(self, indiv):
     """ Removes individual if it exists in history. """
     from pickle import load, dump
     from os.path import exists
     from ..opt import acquire_lock
     if not exists(self.filepath): return 
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: history = load(file)
       if indiv in history: 
         result = history.pop(history.index(indiv))
         with open(self.filepath, "w") as file: dump(history, file)
         return result

   @property 
   def candidates(self):
     """ List of all known candidates. """
     from pickle import load
     from os.path import exists
     from ..opt import acquire_lock
     if not exists(self.filepath): return []
     with acquire_lock(self.filepath) as lock:
       with open(self.filepath, "r") as file: return load(file)

   def __repr__(self):
     """ Python representation of this object. """
     return "{0.__class__.__name__}({1}, {2})".format(self, repr(self.filename), repr(self.directory))
