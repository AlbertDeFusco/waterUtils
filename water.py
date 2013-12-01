#!/usr/bin/env python

import sys
import numpy as np
import math
import molecule

gasReference = np.array(
  [[ 0.0 ,  0.0         ,  0.0650980307 ],
   [ 0.0 ,  0.756950327 , -0.5207842453 ],
   [ 0.0 , -0.756950327 , -0.5207842453 ]] )


def distance(atom1,atom2):
  diff=atom2-atom1
  return np.sqrt((diff**2).sum())

def norm(myArray):
  return np.sqrt((myArray**2).sum())

#probably not very important,
# but this might be made faster
def makeWaters(atomTypes):
  waters=list()
  for oxygen in atomTypes['O']:
    thisWater=water([('O',oxygen)])
    for hydrogen in atomTypes['H']:
      if distance(oxygen,hydrogen) < 1.0:
        thisWater.addAtom(('H',hydrogen))
    thisWater.asArray() #should not be necessary
    waters.append(thisWater)
  return waters

class water(molecule.Molecule):
  def __init__(self,atoms):
    super(water,self).__init__(atoms)

    self.phi=0.0
    self.theta=0.0
    self.psi=0.0

  def getEuler(self):
    self.com = (self.coords[1,:] + self.coords[2,:] + self.coords[0,:]*16.0)/18.0

    rp1 = self.coords[1,:] - self.coords[0,:]
    rp3 = self.coords[2,:] - self.coords[0,:]

    zaxis = rp1 + rp3
    zaxis = -zaxis/norm(zaxis)

    rp1p = sum(rp1*zaxis)

    yaxis = rp1 - rp1p*zaxis
    yaxis = yaxis/norm(yaxis)

    xaxis = np.cross(yaxis,zaxis)

    ## step 2: convert to euler angles
    self.theta = math.acos(zaxis[2])
    sintheta = math.sin(self.theta)

    self.phi=0.0
    self.psi=0.0
    if(abs(sintheta) > 1.e-5):
      tempcos  = -xaxis[2]/sintheta
      acostemp = math.acos(tempcos)
      if(tempcos <= -1.):
        acostemp = math.pi
      if(yaxis[2] >= 0):
        self.psi = acostemp
      else:
        self.psi = 2.0*math.pi - acostemp

      #rare case when zaxis[2] is nearly 1.0
      tempcos = zaxis[0]/sintheta
      tempcos = min(tempcos,1.)
      tempcos = max(tempcos,-1.)

      if(zaxis[1] >=0 ):
        self.phi = math.acos(tempcos)
      else:
        self.phi = 2.0*math.pi - math.acos(tempcos)

    else:
      self.phi=0.0
      if(xaxis[1] >= 0):
        self.psi = math.acos(yaxis[1])
      else:
        self.psi = 2.0*math.pi - math.acos(yaxis[1])

  def immobilize(self,reference):
    #rotate the reference into the phi, theta and psi angles
    c1=math.cos(self.phi)
    c2=math.cos(self.theta)
    c3=math.cos(self.psi)
    s1=math.sin(self.phi)
    s2=math.sin(self.theta)
    s3=math.sin(self.psi)

    # http://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    # the z1y2z3 matrix
    rotate = np.array(
        [ [  c1*c2*c3-s1*s3 , -c3*s1-c1*c2*s3 , c1*s2 ],
          [  c1*s3+c2*c3*s1 ,  c1*c3-s1*c2*s3 , s1*s2 ],
          [ -c3*s2          ,  s2*s3          , c2    ] ] )

    rOxy=reference[0,:].dot(rotate.T) + self.com
    rHy1=reference[1,:].dot(rotate.T) + self.com
    rHy2=reference[2,:].dot(rotate.T) + self.com

    return [('O',rOxy),('H',rHy1),('H',rHy2)]



if __name__ == "__main__":
  #### main #####

  try:
    inputXYZ = open(sys.argv[1],'r')
  except:
    print "Please specify a file name"
    sys.exit(1)

  try:
    nAtoms = int(inputXYZ.readline())
  except ValueError:
    print "The first line is not an integer"
    sys.exit(1)

  comment = inputXYZ.readline()
  #print comment

  #the input geometry is a list of tuples
  inputAtoms = [tuple(thisAtom.split()[0:4]) for thisAtom in inputXYZ]
  inputXYZ.close()

  ## Step 1: dictionary of atoms and coordinates
  atomTypes = molecule.collectAtoms(inputAtoms)

  ## Step 2: assign hydrogens to oxygens
  # The default cutoff is 1 angstrom
  # any other atom type is ignored
  waters=makeWaters(atomTypes)

  ######## Immobilize
  ## step 1: get water-centered axes
  for water in waters:

    water.getEuler()
    water.atoms = water.immobilize(gasReference)
    water.asArray()


  molecule.printCoords(waters,comment)
