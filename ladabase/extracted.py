""" Pre-extracted database subpackage. """

class Encode(object):
  """ Functor for encoding objects. """
  def __init__(self):
    """ Initializes decoding class. """
    from inspect import ismethoddescriptor, isdatadescriptor, isgetsetdescriptor
    from ..vasp.extract import ExtractCommon, ExtractDFT, ExtractGW
    super(Encode, self).__init__()

    def is_descriptor(name):
      return (ismethoddescriptor(name) or isdatadescriptor(name) or isgetsetdescriptor(name))

    self.from_json = {}, {}, {}
    """ Dictionary of values which need be decoded from json. """
    self.sections = {}, {}, {}
    """ Attributes to encode, with section names. """
    for extractor, items, section in zip([ExtractCommon, ExtractDFT, ExtractGW], self.from_json, self.sections):
      for key in dir(extractor):
        if key[0] == '_': continue
        if key == 'comm': continue
        if key == 'directory': continue
        value = getattr(extractor, key)
        if not is_descriptor(value): continue
        value = getattr(value, 'fget', None)
        if value == None: continue
        if not hasattr(value, 'section'): continue
        section[key] = getattr(value, 'section', None)
        json = getattr(value, 'to_json', None)
        if json == None: continue
        items[key] = json 
  
  def __call__(self, extractor, items=None):
    """ Returns dictionary of encoded values. """
    jsons, sections = self.from_json[0].copy(), self.sections[0].copy()
    if extractor.is_dft:
      jsons.update(self.from_json[1])
      sections.update(self.sections[1])
    elif extractor.is_gw:
      jsons.update(self.from_json[2])
      sections.update(self.sections[2])
    if items == None: items = sections.keys()

    result = {}
    for key in items:
      if key not in sections: raise KeyError("Unknown property {0}.".format(key))
      # extracts value.
      try: value = getattr(extractor, key, None)
      except:
        if key == "structure": raise 
        continue
      if value is None: continue

      # stores result.
      if sections[key] == None: section = result
      else:
        if sections[key] not in result: result[sections[key]] = {}
        section = result[sections[key]]
      if key in section:
        raise RuntimeError( "Found two data object with same name {0}{1}."\
                            .format( key, "" if section is result\
                                     else " in section {0}".format(sections[key] )))
      section[key] = jsons[key](value) if key in jsons else value
    
    # add raw data.
    result['raw'] = extractor._id

    return result

class Decode(dict):
  """ Decoding class for pre-extracted vasp-values. """
  def __init__(self):
    """ Initializes decoding class. """
    from inspect import ismethoddescriptor, isdatadescriptor, isgetsetdescriptor
    from ..vasp.extract import ExtractCommon, ExtractDFT, ExtractGW
    super(Decode, self).__init__()

    def is_descriptor(name):
      return (ismethoddescriptor(name) or isdatadescriptor(name) or isgetsetdescriptor(name))

    self.from_json, self.sections = {}, {}, set()
    """ Dictionary of values which need be decoded from json. """
    for extractor in [ExtractCommon, ExtractDFT, ExtractGW]:
      for key in dir(extractor):
        if key[0] == '_': continue
        if key == 'comm': continue
        if key == 'directory': continue
        value = getattr(extractor, key)
        if not is_descriptor(value): continue
        value = getattr(value, 'fget', None)
        if value == None: continue
        section = getattr(value, "section", "Not a section")
        if section == "Not a section": continue
        if section != None: self.sections.add(section)
        value = getattr(value, 'from_json', None)
        if value == None: continue
        self.from_json[key] = value

  def __setitem__(self, key, value):
    """ Sets dictionary items. """
    from .vasp import VaspExtract
    # value has already been set.
    assert key not in self, RuntimeError("Values are read-only.")

    # transform value if needed.
    if key == "raw": value = VaspExtract(value)
    elif key in self.sections:
      for k, v in value.iteritems():
        if key in self.from_json:
          try: v = self.from_json[key](value)
          except: print key, key in self.from_json
        super(Decode, self).__setitem__(k, v)
      return
    elif key in self.from_json:
      try: value = self.from_json[key](value)
      except: print key, key in self.from_json
    super(Decode, self).__setitem__(key, value)

  def __getattr__(self, key):
    """ Links dictionary objects as attributes. """
    assert key in self, AttributeError("Unknown attribute {0}.".format(key))
    return self[key]

  def __setattr__(self, key, value):
    """ Makes dictionary attributes read-only. """
    if key in self: raise AttributeError("Database values are read-only.")
    super(Decode, self).__setattr__(key, value)

  def __dir__(self):
    """ Returns list of attributes. """
    return list(   set([k for k in dir(self.__class__) if k[0] != '_']) \
                 | set([k for k in self.keys() if k[0] != '_']) \
                 | set(['_id']) )

def to_secondary(collection='extracted', filter=None, items=None, update=False):
  """ Extracts value to secondary database. """
  from . import Manager
  from .vasp import VaspExtract
  ladabase = Manager()
  encoder = Encode()
  collection = ladabase.database[collection]
  for element in ladabase.files.find(filter):
    extract = VaspExtract(element)
    encoded = encoder(extract)

    indatabase = [k for k in collection.find({'_id': extract._id})]
    if len(indatabase) == 1:
      if update: encoded['_id'] = indatabase[0]['_id']
      else: continue
    elif len(indatabase) > 1:
      raise RuntimeError('found more than one items corresponding to the same object {0}.'.format(extract._id))
    collection.save(encoded)
