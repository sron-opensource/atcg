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
GRID,5,2,50.,0.,-50.,3
GRID,6,2,50.,0.,50.,3
CQUAD4,1,1,6,5,1,2,,0.,
CQUAD4,2,1,5,6,4,3,,0.,
$ ----------------------------------------
SPC,161,1,,2.
SPC,161,2,,2.
SPC,161,3,,1.
SPC,161,4,,1.
ENDDATA
