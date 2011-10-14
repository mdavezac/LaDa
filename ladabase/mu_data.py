""" Elemental mus (FERE). """
from quantities import eV

elemental_mus = { 'Ag': -0.82700958541595615*eV,
                  'Au': -2.2303410086960551*eV,
                  'Ba': -1.3944992462870172*eV,
                  'Be': -3.3972092621754264*eV,
                  'Bi': -4.3853003286558812*eV,
                  'Cl': -1.6262437135301639*eV,
                  'Co': -4.7543486260270402*eV,
                  'Cr': -7.2224146752384204*eV,
                  'Cu': -1.9725806522979044*eV,
                  'F' : -1.7037867766570287*eV,
                  'Fe': -6.1521343161090325*eV,
                  'Ge': -4.137439286830797*eV,
                  'Hf': -7.397695761161847*eV,
                  'Hg': -0.12361566177444684*eV,
                  'Ir': -5.964577394407752*eV,
                  'K' : -0.80499202755075006*eV,
                  'La': -3.6642174822805287*eV,
                  'Li': -1.6529591953887741*eV,
                  'Mn': -6.9965778258511993*eV,
                  'Na': -1.0640326227725869*eV,
                  'Nb': -6.6867516375690608*eV,
                  'Ni': -3.5687859474688026*eV,
                  'Pd': -3.1174044624888873*eV,
                  'Pt': -3.9527597082085424*eV,
                  'Rb': -0.6750560483522855*eV,
                  'Rh': -4.7622899695820369*eV,
                  'Sb': -4.2862260747305099*eV,
                  'Sc': -4.6302422200922519*eV,
                  'Si': -4.9927748122726356*eV,
                  'Sn': -3.7894939351245469*eV,
                  'Sr': -1.1674559193419329*eV,
                  'Ta': -8.8184831379805324*eV,
                  'Te': -3.2503408197224912*eV,
                  'Ti': -5.5167842601434147*eV,
                  'V' : -6.4219725884764864*eV,
                  'Y' : -4.812621315561298*eV,
                  'Zr': -5.8747056261113126*eV,
                  'O' :-4.76*eV,
                  'S' :-4.00*eV,
                  'Se':-3.55*eV,
                  'N' :-8.51*eV,
                  'P' :-5.64*eV,
                  'As':-5.06*eV,
                  'Mg':-0.99*eV,
                  'Ca':-1.64*eV,
                  'Zn':-0.84*eV,
                  'Cd':-0.56*eV,
                  'Ga':-2.37*eV,
                  'Al':-3.02*eV,
                  'In':-2.31*eV}
""" mus determined from original and FERE paper. """


original_paper = ['O', 'S', 'Se', 'N', 'P', 'As', 'Mg', 'Ca', 'Zn', 'Cd', 'Ga', 'Al', 'In']
""" Atoms determined from original Stephan Lany paper. """

def enthalpy(extract):
  """ Returns enthalpy if compuational parameters are ok, else None. """
  from quantities import eV
  from .. import periodic_table as pt

  if not extract.is_dft: return False, False
  if extract.pseudotential != 'PAW_PBE': return False, False
  if extract.encut.rescale(eV).magnitude > 339.0: return False, False
  for specie in extract.species:
    if specie not in elemental_mus: return None
    U = None
    if specie in ["Cu", "Ag"]: U = 5 # in eV
    elif pt.__dict__[specie].row >= 3 or pt.__dict__[specie].row <= 11: U = 3 # in eV
    if U == None and specie in extract.HubbardU_NLEP: return False, False
    elif U != None:
      if specie not in extract.HubbardU_NLEP: return False, False
      if len(extract.HubbardU_NLEP[specie]) != 1: return False, False
      if extract.HubbardU_NLEP["type"] != "dudarev": return False, False
      if extract.HubbardU_NLEP["l"] != 2: return False, False
      if extract.HubbardU_NLEP["func"] != "U": return False, False
      if abs(extract.HubbardU_NLEP["J"]) > 1e-3: return False, False
      if abs(extract.HubbardU_NLEP["U"]-U) > 1e-3: return False, False
  result = extract.energy - sum([n * elemental_mus[s] for n, s in zip(extract.stoechiometry, extract.species)])
  return result / float(sum(extract.stoechiometry))

