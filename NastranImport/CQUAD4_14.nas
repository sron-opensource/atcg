$ Finite Element Mesh - PTC Inc. - MSC/NASTRAN 2012 - 000-017-1017
$
SOL SESTATIC
CEND
DISPLACEMENT = ALL
SPCFORCES = ALL
FORCE = ALL
OLOAD = ALL
STRESS(SORT1,REAL,VONMISES,BILIN) = ALL
STRAIN(SORT1,REAL,VONMISES,BILIN,FIBER) = ALL
ESE = ALL
BEGIN BULK
PARAM,POST,0
PARAM,AUTOSPC,YES
$ Global Coordinate System of the model
CORD2R,1,0,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
$ ----------------------------------------
$ Mesh "000-017-1017"
$ Included Components :
$   000-017-1017
$ Coordinate System for grid coordinates
CORD2R,2,1,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
$ Coordinate System for default grid displacement
CORD2R,3,2,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
MAT1,1,0.,,0.,1.E-10,0.,,,
PSHELL,1,1,5.,1,,1
GRID,1,2,-25.,5.,0.,3
GRID,2,2,-25.,-5.,0.,3
GRID,3,2,25.,5.,0.,3
GRID,4,2,25.,-5.,0.,3
GRID,5,2,-25.,0.,0.,3
GRID,6,2,20.4545,5.,0.,3
GRID,7,2,13.6364,5.,0.,3
GRID,8,2,4.54545,5.,0.,3
GRID,9,2,-4.54545,5.,0.,3
GRID,10,2,-13.6364,5.,0.,3
GRID,11,2,-20.4545,5.,0.,3
GRID,12,2,25.,0.,0.,3
GRID,13,2,-20.4545,-5.,0.,3
GRID,14,2,-13.6364,-5.,0.,3
GRID,15,2,-4.54545,-5.,0.,3
GRID,16,2,4.54545,-5.,0.,3
GRID,17,2,13.6364,-5.,0.,3
GRID,18,2,20.4545,-5.,0.,3
GRID,19,2,-20.4545,0.,0.,3
GRID,20,2,-13.6364,0.,0.,3
GRID,21,2,-4.54545,0.,0.,3
GRID,22,2,4.54545,0.,0.,3
GRID,23,2,13.6364,0.,0.,3
GRID,24,2,20.4545,0.,0.,3
CQUAD4,1,1,24,12,3,6,,0.,
CQUAD4,2,1,23,24,6,7,,0.,
CQUAD4,3,1,22,23,7,8,,0.,
CQUAD4,4,1,21,22,8,9,,0.,
CQUAD4,5,1,20,21,9,10,,0.,
CQUAD4,6,1,19,20,10,11,,0.,
CQUAD4,7,1,5,19,11,1,,0.,
CQUAD4,8,1,18,4,12,24,,0.,
CQUAD4,9,1,17,18,24,23,,0.,
CQUAD4,10,1,16,17,23,22,,0.,
CQUAD4,11,1,15,16,22,21,,0.,
CQUAD4,12,1,14,15,21,20,,0.,
CQUAD4,13,1,13,14,20,19,,0.,
CQUAD4,14,1,2,13,19,5,,0.,
$ ----------------------------------------
ENDDATA
