import openbabel
import sys
import pickle
#import plpy
def notice(msg):
 print msg
 plpy.notice(msg)

# functions used in atoms, bonds and fragments functions
tbl = openbabel.OBElementTable()

def nbr_smarts(a1, a2=None):
 # construct smarts of each neighbor, possible excluding a2
 smarts = [] # each neighbor as array element for later sort
 for nbr in openbabel.OBAtomAtomIter(a1):
  if nbr.GetAtomicNum() == 1: continue
  #if a2 and nbr.GetIdx() == a2.GetIdx():
  # continue
  # include bond order if necessary (not single or aromatic)
  b = a1.GetBond(nbr)
  if b.IsSingle() and a1.IsAromatic():
   bnd = "-"
  elif b.IsDouble():
   bnd = "="
  elif b.IsTriple():
   bnd = "#"
  else:
   bnd = ""
  smarts.append(bnd + atom_sym(nbr,True))
 # now sort each neighbor smarts to provide unique ordering
 smarts.sort()
 slen = len(smarts)
 sma = ""
 # assemble final smarts, if there are neighbors
 if slen > 0:
 # include parens when 2 or more neighbors
  for i in range(0, slen-1):
   sma += "(" + smarts[i] + ")"
  sma += smarts[slen-1]
 return sma

def atom_sym(a, smarts_flag=False):
 "construct atomic smarts with symbol, charge, hcount, connect count"
 atnum = a.GetAtomicNum()
 if atnum > 0:
  atsym = tbl.GetSymbol(atnum)
 else:
  atsym = "*"
 if a.IsAromatic(): atsym = atsym.lower()
 if smarts_flag:
  s = "[" + atsym
  s += "%+d" % a.GetFormalCharge()
  s += "h%d" % a.ImplicitHydrogenCount()
  s += "D%d" % (a.GetHvyValence() + a.ExplicitHydrogenCount())
  s += "]"
 else:
  if atnum in [1,5,6,7,8,9,15,16,17,35,53]:
   s = atsym
  else:
   s = '[' + atsym + ']'
 return s

def bond_sym(a,b):
 s = ""
 bond = a.GetBond(b)
 if bond:
  if   bond.IsDouble(): s += "="
  elif bond.IsTriple(): s += "#"
 return s

def apath(alist, smarts_flag):
#s = "".join([atom_sym(a) for a in alist])
 """return a string of atom symbols from atoms in alist,
    with appropriate bond symbols
 """
 a = alist[0]
 s = atom_sym(a, smarts_flag)
 for i in range(1, len(alist)):
  b = alist[i]
  s += bond_sym(a,b) + atom_sym(b, smarts_flag)
  a = b
 return s
 
def path(alist, smarts_flag):
 """return string (smarts) representation of atoms in alist
    special case for rings.
 """
 # check for a ring
 if len(alist) > 2 and alist[0].IsConnected(alist[-1]):
   rpath = ring_path(alist, smarts_flag)
 else:
   rpath = None

 # make path and its reverse.  return lesser of the two
 s = apath(alist, smarts_flag)
 alist.reverse()
 r = apath(alist, smarts_flag)
 if s < r: return s,rpath
 else:     return r,rpath

def aring_path(alist, start, smarts_flag):
 """string representation of a ring of atoms in alist
  starting with atom #start, going to end, then looping
  back to 0 and up to start"""
 a = alist[start]
 s = atom_sym(a, smarts_flag) + '1'
 for i in range(start+1, len(alist)):
  b = alist[i]
  s += bond_sym(a,b) + atom_sym(b, smarts_flag)
  a = b
 for i in range(0, start):
  b = alist[i]
  s += bond_sym(a,b) + atom_sym(b, smarts_flag)
  a = b
 return s + '1'
  
def ring_path(alist, smarts_flag):
 """return string representation of ring of atoms in alist
    taking care to examine all ways ring can be formed, starting
    with each atom or reversed"""
 s = list()
 for j in range(0, len(alist)):
  s.append(aring_path(alist, j, smarts_flag))
 alist.reverse()
 for j in range(0, len(alist)):
  s.append(aring_path(alist, j, smarts_flag))
 s.sort()
 return s[0]

def inlist(atom, alist):
 return [a.GetIdx() for a in alist].count(atom.GetIdx()) > 0
#ndx = atom.GetIdx()
#d = [a.GetIdx() for a in alist]
#return d.count(ndx) > 0

def make_atom_paths(molpath, alist, from_atom, atom, depth, minpath, maxpath, smarts_flag):
 if depth >= maxpath:
  return depth

 alist.append(atom)
 for nbr in openbabel.OBAtomAtomIter(atom):
  if nbr.GetIdx() == from_atom.GetIdx(): continue # don't go backwards where we came from
  if nbr.GetAtomicNum() < 2: continue # break path for * or H atom
  if inlist(nbr, alist):
   #alist.append(nbr)
   break
  else:
   nbrlist = list()
   nbrlist.extend(alist)
   make_atom_paths(molpath, nbrlist, atom, nbr, depth+1, minpath, maxpath, smarts_flag)

 if depth+1 >= minpath:
  spath,rpath = path(alist, smarts_flag)
  if spath: molpath.append(spath)
  if rpath: molpath.append(rpath)

def make_mol_paths(mol, minpath, maxpath, smarts_flag):
 alist = list()
 molpath = list()
 for a in openbabel.OBMolAtomIter(mol):
  if a.GetAtomicNum() < 2: continue # break path for * or H atom
  make_atom_paths(molpath,alist,a,a,0,minpath,maxpath,smarts_flag)
  del alist[:]
 return molpath

def mol2amol(mol):
  mol.Kekulize()
  atoms = list()
  for a in openbabel.OBMolAtomIter(mol):
    atom = list()
    atom.append(a.GetHyb())
    atom.append(a.GetAtomicNum())
    atom.append(a.GetFormalCharge())
    atom.append(a.GetIsotope())
    if a.IsClockwise():
      stereo = 1
    elif a.IsAntiClockwise():
      stereo = 2
    elif a.IsChiral():
      stereo = 3
    else:
      stereo = 0
    atom.append(stereo)
    atom.append(a.GetSpinMultiplicity())
    if a.IsAromatic():
      aromatic = 1
    else:
      aromatic = 0
    atom.append(aromatic)

    atoms.append(atom)

  bonds = list()
  for b in openbabel.OBMolBondIter(mol):
    bond = list()
    bond.append(b.GetBeginAtomIdx())
    bond.append(b.GetEndAtomIdx())
    bond.append(b.GetBondOrder())
    if b.IsWedge():
      stereo = 1
    elif b.IsHash():
      stereo = 6
    else:
      stereo = 0
    bond.append(stereo)
    if b.IsAromatic():
      aromatic = 1
    else:
      aromatic = 0
    bond.append(aromatic)
    if b.IsUp():
      updown = 1
      print updown
    elif b.IsDown():
      updown = 2
      print updown
    else:
      updown = 0
    bond.append(updown)

    bonds.append(bond)

  return pickle.dumps([atoms, bonds],2)
  #return [atoms, bonds]
 

def amol2mol(amol):
  mol = openbabel.OBMol()
  atoms = list()
  bonds = list()
  atoms,bonds = pickle.loads(amol)
  natoms = len(atoms)
  mol.ReserveAtoms(natoms)
  atomidx = 0
  chiatoms = dict()
  a = openbabel.OBAtom()
  for atom in atoms:
    hybridization, atnum, fcharge, isotope, stereo, mult, aromatic = atom
    a.Clear()
    atomidx += 1
    a.SetIdx(atomidx)
    a.SetHyb(hybridization)
    a.SetAtomicNum(atnum)
    a.SetIsotope(isotope)
    a.SetFormalCharge(fcharge)
    if stereo == 1:
      a.SetClockwiseStereo()
    elif stereo == 2:
      a.SetAntiClockwiseStereo()
    elif stereo == 3:
      a.SetChiral()
    a.SetSpinMultiplicity(mult)
    if aromatic: a.SetAromatic()

    if not mol.AddAtom(a):
      return None

    if stereo:
      cd = openbabel.OBChiralData()
      catom = mol.GetAtom(atomidx)
      catom.CloneData(cd)
      chiatoms[catom] = cd
      pass

  for bond in bonds:
    flags = 0
    #b = openbabel.OBBond()
    start, end, order, stereo, aromatic, updown = bond
    if start==0 or end==0 or order==0 or start>natoms or end>natoms:
      return None
    if order == 4: order = 5
    if stereo == 1:
      flags |= openbabel.OB_WEDGE_BOND
    elif stereo == 6:
      flags |= openbabel.OB_HASH_BOND
    if aromatic:
      flags |= openbabel.OB_AROMATIC_BOND
    if updown == 1:
      flags |= openbabel.OB_TORUP_BOND
      print flags
    elif updown == 2:
      flags |= openbabel.OB_TORDOWN_BOND
      print flags
    if not mol.AddBond(start, end, order, flags):
      return None

    chidata = chiatoms.get(mol.GetAtom(start))
    if chidata:
      chidata.AddAtomRef(end, openbabel.input)
    chidata = chiatoms.get(mol.GetAtom(end))
    if chidata:
      chidata.AddAtomRef(start, openbabel.input)

  mol.SetAromaticPerceived()
  mol.SetKekulePerceived()
  return mol
  
class obchord:

  def __init__(self, Global):
    """use Global dictionary supplied by plpython to buffer parsed
       smiles and smarts for quicker lookup when reused.
    """
    # buffering smiles and smarts can eat up lots of memory,
    #  but speeds things up considerably
    #  certainly if an entire table can acutally fit
    # but also if this is the maximum size of a hit list, say from match
    if not Global.has_key("OB"): Global["OB"] = dict()
    self.GOB = Global["OB"]
    self.maxsmi = 1000
    self.maxqsmi = 1000
    if not self.GOB.has_key("mol"): self.GOB["mol"] = dict()
    self.mol = self.GOB["mol"]
    # pick a reasonable number of smarts patterns you expect to use often,
    # say 166 public keys or even 1000 fragment keys
    self.maxsma = 1000
    if not self.GOB.has_key("pat"): self.GOB["pat"] = dict()
    self.pat = self.GOB["pat"]
    if not self.GOB.has_key("obc"): self.GOB["obc"] = openbabel.OBConversion()
    self.obc = self.GOB["obc"]
    self.obc.SetInFormat("smi")
#   if not self.GOB.has_key("fobc"): self.GOB["fobc"] = openbabel.OBConversion()
#   self.fobc = self.GOB["fobc"]
#   self.fobc.SetInAndOutFormats("sdf","can")
    if not self.GOB.has_key("fmol"): self.GOB["fmol"] = openbabel.OBMol()
    self.fmol = self.GOB["fmol"]
    if not self.GOB.has_key('fp'): self.GOB['fp'] = openbabel.vectorUnsignedInt()
    if not self.GOB.has_key('fpr'): self.GOB['fpr'] = openbabel.OBFingerprint.FindFingerprint('FP2')
    self.fp = self.GOB["fp"]
    self.fpr = self.GOB["fpr"]
    if not self.GOB.has_key("isomapper"): self.GOB["isomapper"] = dict()
    self.isomapper = self.GOB["isomapper"]
    self.isomap = openbabel.vpairUIntUInt()

  def smilesBuffer(self,smi):
    """one, or all keys(smiles) stored in global
    """
    if smi:
      return [self.mol[smi]]
    else:
      return self.mol.keys()
  
  def smartsBuffer(self,sma):
    """one, or all keys(smarts) stored in global
    """
    if sma:
      return [self.pat[sma]]
    else:
      return self.pat.keys()
  
  def parse_smi(self,smi):
    """parse smiles and return OBMol after storing in global dict
       or return from global dict"""
    if self.mol.has_key(smi):
      # return copy is slower, but safer?
      #return openbabel.OBMol(self.mol[smi])
      #notice('found mol for %s' % smi)
      return self.mol[smi]
    if len(self.mol) < self.maxsmi:
      amol = openbabel.OBMol()
      #notice('new mol for %s' % smi)
    else:
      key,amol = self.mol.popitem()
      #notice('mol reuse %s for %s' % (key,smi))
      amol.Clear()
    if self.obc.ReadString(amol, smi):
      self.mol[smi] = amol
      # return copy is slower, but safer?
      # return openbabel.OBMol(amol)
      return amol
    else:
      return None

  def parse_qsmi(self,qsmi):
    """parse query smiles and return OBIsomorphismMapper instance after storing in global dict
       or return from global dict"""
    if self.isomapper.has_key(qsmi):
      return self.isomapper[qsmi]
    query = openbabel.CompileSmilesQuery(qsmi)
    if query:
      if len(self.isomapper) < self.maxqsmi:
        #notice('new query for %s' % qsmi)
        pass
      else:
        key,mapper = self.isomapper.popitem()
        #notice('query reuse %s for %s' % (key,qsmi))
        del mapper
      mapper = openbabel.OBIsomorphismMapper.GetInstance(query)
      self.isomapper[qsmi] = mapper
      return mapper
    else:
      return None
  
  def parse_sma(self,sma):
    """parse smarts and return OBSmartsPattern after storing in global dict
       or return from global dict"""
    if self.pat.has_key(sma):
      return self.pat[sma]
    if len(self.pat) < self.maxsma:
      pat = openbabel.OBSmartsPattern()
      #notice('new pat for %s' % sma)
    else:
      key,pat = self.pat.popitem()
      #notice('pattern reuse %s for %s' % (key,sma))
    if pat.Init(sma):
      self.pat[sma] = pat
      return pat
    else:
      #notice('pattern None')
      return None
  
  def input_formats(self):
    return self.obc.GetSupportedInputFormat()

  def output_formats(self):
    return self.obc.GetSupportedOutputFormat()

  def parse_format(self,strfile,fmt):
    """parse string formatted as fmt and return mol
    """
    fobc = openbabel.OBConversion()
    if fobc.SetInFormat(fmt):
      try:
        self.fmol.Clear();
        if fobc.ReadString(self.fmol, strfile):
          return self.fmol
      except:
        raise ValueError("Error parsing input file")
    else:
      raise ValueError("Unrecognized input format type")
    return None
  
  def writestring(self,mol,fmt,opt=None):
    """convert mol to output format and return string"""
    fobc = openbabel.OBConversion()
    if fobc.SetOutFormat(fmt):
      try:
        #if opt: fobc.AddOption(opt,fobc.OUTOPTIONS)
        if opt: fobc.SetOptions(opt,fobc.OUTOPTIONS)
        outstr = fobc.WriteString(mol,1)
        #notice('length of output %s' % len(outstr))
        #if opt: fobc.RemoveOption(opt,fobc.OUTOPTIONS)
        return outstr
      except:
        raise ValueError("Error writing ouput file")
    else:
      raise ValueError("Unrecognized output format type")
    return None

  def fingerprint(self,mol,nbits):
    # using static in Global saves lots of allocation/free execution time
    for i in range(0,len(self.fp)):self.fp[i] = 0
    self.fpr.GetFingerprint(mol,self.fp,nbits)
    return "".join(["%08x" % x for x in self.fp])

  def substructure(self, amol, atoms):
    # return substructure of mol containing only atoms
    delatom = list()
    # copy amol to alter it
    omol = openbabel.OBMol(amol)
    for idx in range(omol.NumAtoms()):
      if not idx in atoms: delatom.append(omol.GetAtom(1+idx))
    for a in delatom:
      omol.DeleteAtom(a)
    return omol

  def contains(self,mol,mapper):
    # returns true is mol is a superstructure of/contains qsmi
    # using static in Global saves lots of allocation/free execution time
    #isomap = openbabel.vpairUIntUInt()
    del self.isomap[0:]
    mapper.MapFirst(mol,self.isomap)
    if self.isomap: return True
    else: return False
