$ Finite Element Mesh - PTC Inc. - MSC/NASTRAN 2012 - 000-016-1364
$
SOL SESTATIC
APP HEAT
CEND
ANALYSIS = HEAT
THERMAL = ALL
FLUX = ALL
SPCFORCES = ALL
OLOAD = ALL
SUBCASE = 1
   SPC = 161
BEGIN BULK
PARAM,POST,0
PARAM,AUTOSPC,YES
$ Global Coordinate System of the model
CORD2R,1,0,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
$ ----------------------------------------
$ Mesh "000-016-1364"
$ Included Components :
$   000-016-1364
$ Coordinate System for grid coordinates
CORD2R,2,1,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
$ Coordinate System for default grid displacement
CORD2R,3,2,0.,0.,0.,0.,0.,1.,
,1.,0.,0.
MAT4,1,236.,8.93E8,2.7E-9
PSHELL,1,1,3.,1,,1
GRID,1,2,-50.,0.,-50.,3
GRID,2,2,-50.,0.,50.,3
GRID,3,2,150.,0.,-50.,3
GRID,4,2,150.,0.,50.,3
GRID,5,2,110.,0.,10.,3
GRID,6,2,80.,0.,-20.,3
GRID,7,2,-5.,0.,5.,3
GRID,8,2,30.,0.,25.,3
GRID,9,2,-50.,0.,0.,3
GRID,10,2,100.,0.,-50.,3
GRID,11,2,50.,0.,-50.,3
GRID,12,2,0.,0.,-50.,3
GRID,13,2,150.,0.,10.,3
GRID,14,2,-5.,0.,50.,3
GRID,15,2,50.,0.,50.,3
GRID,16,2,100.,0.,50.,3
GRID,17,2,110.,0.,-19.4351,3
GRID,18,2,30.,0.,50.,3
GRID,19,2,67.5,0.,16.25,3
GRID,20,2,150.,0.,-17.7405,3
CQUAD4,1,1,8,18,15,19,,0.,
CQUAD4,2,1,6,11,8,19,,0.,
CQUAD4,3,1,5,19,15,16,,0.,
CQUAD4,4,1,17,20,3,10,,0.,
CQUAD4,5,1,6,17,10,11,,0.,
CQUAD4,6,1,5,13,20,17,,0.,
CQUAD4,7,1,5,16,4,13,,0.,
CQUAD4,8,1,5,17,6,19,,0.,
CQUAD4,9,1,7,9,2,14,,0.,
CQUAD4,10,1,14,18,8,7,,0.,
CQUAD4,11,1,9,7,12,1,,0.,
CQUAD4,12,1,7,8,11,12,,0.,
$ ----------------------------------------
SPC,161,1,,2.
SPC,161,2,,2.
SPC,161,3,,1.
SPC,161,4,,1.
SPC,161,9,,2.
SPC,161,13,,1.
SPC,161,20,,1.
ENDDATA
