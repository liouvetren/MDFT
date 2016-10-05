source leaprc.ff03.r1
source leaprc.gaff
set default pbradii mbondi
lig = loadmol2 lig.ac.mol2
frcmod = loadamberparams lig.leap.frcmod
saveamberparm lig lig.leap.prm lig.leap.crd
savepdb lig lig.leap.pdb
quit
