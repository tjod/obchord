import sys
import obchord
import openbabel

def obmol(imol, mol):
  if mol:
    name = mol.GetTitle()
    mol.SetTitle("")
    if name == "": name = "\\N"
    print str(imol) + "\t" + name + "\t" + oc.writestring(mol,"can","i") + "\t" + oc.writestring(mol,"can").replace("\\","\\\\") + "\t\\X" + oc.fingerprint(mol, 1024)
  else:
    print str(imol) + "\t\\N\t\\N\t\\N\t\\N"

def obprop(imol, mol):
 for p in mol.GetData():
  if openbabel.toPairData(p).GetDataType() == openbabel.PairData:
    if p.GetAttribute() == 'OpenBabel Symmetry Classes':
      pass
    else:
      print str(imol) + "\t" + p.GetAttribute() + "\t" +  p.GetValue()

def newCopy(properties):
  if properties:
    print "Copy properties (id, name, value) From Stdin;"
  else:
    print "Copy ob_structure (id, name, cansmiles, isosmiles, fp) From Stdin;"

properties = False
nchunk = 1000
if len(sys.argv) > 2:
  schema = sys.argv[1];
  molfile = sys.argv[2]
  if len(sys.argv) > 3:
    if sys.argv[3] == "--properties":
      properties = True
      nchunk = 10000
  else:
    properties = False
    nchunk = 1000
else:
  print >> sys.stderr, "usage: python sdfloader.py schema sdfile [--properties]\n; writes to stdout intended for psql input"
  exit(0)

print """\\timing
--Drop Schema If Exists %s Cascade;
Create Schema %s;
Grant Usage On Schema %s To Public;
Set search_path=%s;""" % (schema, schema, schema, schema)

print """
Create Table ob_structure (
 id Integer Primary Key,
 name Text,
 cansmiles Text,
 isosmiles Text,
 fp Bit Varying,
 fpbits Int
 );
Create Table properties (
 id Integer,
 name Text,
 value Text);
"""

newCopy(properties)

GD = dict()
oc = obchord.obchord(GD)
molblock = list()
imol = 0
obconversion = openbabel.OBConversion()
obconversion.SetInFormat("sdf")
mol = openbabel.OBMol()
notatend = obconversion.ReadFile(mol, molfile)
while notatend:
    mol.Clear()
    notatend = obconversion.Read(mol)
    imol += 1
    if properties:
      obprop(imol, mol)
    else:
      obmol(imol, mol)
    molblock = list()
    if 0 == (imol % nchunk):
      print >> sys.stderr, imol,
      print "\\."
      newCopy(properties)
