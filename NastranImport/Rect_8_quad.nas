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
GRID,5,2,-50.,0.,0.,3
GRID,6,2,100.,0.,-50.,3
GRID,7,2,50.,0.,-50.,3
GRID,8,2,0.,0.,-50.,3
GRID,9,2,150.,0.,0.,3
GRID,10,2,0.,0.,50.,3
GRID,11,2,50.,0.,50.,3
GRID,12,2,100.,0.,50.,3
GRID,13,2,0.,0.,0.,3
GRID,14,2,50.,0.,0.,3
GRID,15,2,100.,0.,0.,3
CQUAD4,1,1,15,9,3,6,,0.,
CQUAD4,2,1,14,15,6,7,,0.,
CQUAD4,3,1,13,14,7,8,,0.,
CQUAD4,4,1,5,13,8,1,,0.,
CQUAD4,5,1,12,4,9,15,,0.,
CQUAD4,6,1,11,12,15,14,,0.,
CQUAD4,7,1,10,11,14,13,,0.,
CQUAD4,8,1,2,10,13,5,,0.,
$ ----------------------------------------
SPC,161,1,,2.
SPC,161,2,,2.
SPC,161,3,,1.
SPC,161,4,,1.
SPC,161,5,,2.
SPC,161,9,,1.
ENDDATA
