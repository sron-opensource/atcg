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
GRID,13,2,150.,0.,0.,3
GRID,14,2,0.,0.,50.,3
GRID,15,2,50.,0.,50.,3
GRID,16,2,100.,0.,50.,3
GRID,17,2,118.,0.,-22.,3
CTRIA3,1,1,15,8,14,,0.,,
CTRIA3,2,1,15,6,8,,0.,,
CTRIA3,3,1,15,5,6,,0.,,
CTRIA3,4,1,15,16,5,,0.,,
CTRIA3,5,1,10,17,3,,0.,,
CTRIA3,6,1,10,6,17,,0.,,
CTRIA3,7,1,10,11,6,,0.,,
CTRIA3,8,1,5,13,17,,0.,,
CTRIA3,9,1,5,4,13,,0.,,
CTRIA3,10,1,5,16,4,,0.,,
CTRIA3,11,1,5,17,6,,0.,,
CTRIA3,12,1,14,7,2,,0.,,
CTRIA3,13,1,14,8,7,,0.,,
CTRIA3,14,1,9,2,7,,0.,,
CTRIA3,15,1,9,12,1,,0.,,
CTRIA3,16,1,9,7,12,,0.,,
CTRIA3,17,1,13,3,17,,0.,,
CTRIA3,18,1,8,11,7,,0.,,
CTRIA3,19,1,8,6,11,,0.,,
CTRIA3,20,1,12,7,11,,0.,,
$ ----------------------------------------
SPC,161,1,,2.
SPC,161,2,,2.
SPC,161,3,,1.
SPC,161,4,,1.
SPC,161,9,,2.
SPC,161,13,,1.
ENDDATA
