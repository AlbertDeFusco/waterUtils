#!/usr/bin/env python

import numpy as np
import math

def distance(atom1,atom2):
  diff=atom2-atom1
  return np.sqrt((diff**2).sum())

def makeArray(myTuple):
  return np.array([float(i) for i in myTuple])

def collectAtoms(inputXYZ):
  #inputXYZ is a list of tuples of characters
  # Name,x,y,z
  atomTypes=dict()
  for thisAtom in inputXYZ:
    if thisAtom[0] in atomTypes:
      atomTypes[thisAtom[0]].append(makeArray(thisAtom[1:4]))
    else:
      atomTypes[thisAtom[0]] = [makeArray(thisAtom[1:4])]
  return atomTypes

def printCoords(molecules,comment=''):
  numAtoms=sum([len(molecule.atoms) for molecule in molecules])
  print numAtoms
  print comment.strip()
  for mol in molecules:
    for atom in mol.atoms:
      print "%-2s       %12.7f  %12.7f  %12.7f" % (atom[0],atom[1][0],atom[1][1],atom[1][2])


class Molecule(object):
  def __init__(self,atoms=list()):
    #atoms is a list of tuples (Name,np.array(x,y,z))
    self.atoms=atoms
    self.coords=np.array([atm[1] for atm in atoms])
  def addAtom(self,newAtom):
    #newAtom is a tuple (Name,np.array(x,y,z))
    self.atoms.append(newAtom)
    self.coords=np.vstack((self.coords,newAtom[1]))
  def asArray(self):
    self.coords=np.array([atm[1] for atm in self.atoms])
