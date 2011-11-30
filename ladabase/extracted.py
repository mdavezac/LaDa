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
        section[key] = getattr(value, 'section', None)
        if getattr(value, 'to_json', None) is not None:
          items[key] = getattr(value, 'to_json') 
  
  def __call__(self, extractor, items):
    """ Returns dictionary of encoded values. """
    jsons, sections = self.from_json[0].copy(), self.sections[0].copy()
    if extractor.is_dft:
      jsons.update(self.from_json[1])
      sections.update(self.sections[1])
    elif extractor.is_gw:
      jsons.update(self.from_json[2])
      sections.update(self.sections[2])

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


def to_secondary(collection='extracted', filter=None, items=None, update=False, fromextracted=False):
  """ Extracts value to secondary database. """
  from .. import periodic_table as pt
  from ..crystal.read_write import castep
  from . import Manager
  from .vasp import VaspExtract
  from .mu_data import enthalpy
  ladabase = Manager()
  encoder = Encode()
  collection = ladabase.database[collection]
  iterator = ladabase.files.find(filter) if not fromextracted else collection.find(filter)
  for element in iterator:

    if fromextracted:
      encoded, element = element, ladabase.files.find_one({'_id': element['metadata']['raw']}) 
    else: 
      encoded = collection.find_one({'metadata.raw': element['_id']})
      if (not update) and encoded != None: continue
      elif encoded == None: encoded = {'input': {}, 'output': {}, 'metadata': {}}
      
    extract = VaspExtract(element)
    # create dictionary with computational details.
    encoded['input'] = encoder(extract, ['ispin', 'algo', 'encut',
                                         'nelect', 'HubbardU_NLEP', 'pseudopotential'])['input'] 
    try: value = extract.reciprocal_volume.magtnitude.tolist() / float(sum(extract.multiplicity))
    except: pass
    else: encoded['input']['kpoint_density'] = value
    encoded['input']['corrections'] = extract.HubbardU_NLEP
    # create dictionary with output.
    encoded['output'] = encoder(extract, ['total_energy', 'vbm', 'cbm', 'pressure'])['output']
    encoded['output']['total_energy'] /= float(sum(extract.stoechiometry))
    encoded['output']['gap'] = encoded['output']['cbm'] - encoded['output']['vbm']
    try: value = extract.density.magnitude.tolist()
    except: pass
    else: encoded['output']['density'] = value
    # creates dictionary with metadata.
    encoded['metadata']['raw'] = extract._id
    if not (extract.is_dft or extract.is_gw):  encoded['metadata']['functional'] = 'unknown'
    else: encoded['metadata']['functional'] = 'dft' if extract.is_dft else 'gw'
    value = enthalpy(extract)
    if value is not None: encoded['metadata']['Enthalpy'] = value
    encoded['metadata']['structure'] = castep(extract.structure)

    encoded['metadata']['formula'] = ""
    encoded['input']['species'] = {}
    def sortme(args):
      specie = pt.__dict__[args[0]]
      return (specie.column + specie.row * 0.01) if args[0] in pt.__dict__ else 20
    for a, n in sorted([(a, n) for a, n in zip(extract.species, extract.stoechiometry)], key=sortme):
      encoded['metadata']['formula'] += "{0}{1}".format(a, n)
      encoded['input']['species'][a] = n

    collection.save(encoded)
