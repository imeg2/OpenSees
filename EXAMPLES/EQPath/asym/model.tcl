set systemTime [clock seconds]
puts "Starting Analysis: [clock format $systemTime -format "%d-%b-%Y %H:%M:%S"]"
set startTime [clock clicks -milliseconds];
# SET UP ----------------------------------------------------------------------------
wipe;                             # clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6; # Define the model builder, ndm=#dimension, ndf=#dofs

source LibUnits.tcl;                      # define units
set length 1400
set nl 40
set dl [expr $length/$nl]

set StartNode 1;
set EndNode [expr $nl];

#Nodes, NodeNumber, xCoord, yCoord, zCoord
for {set i 1} {$i<=$nl} {incr i 1} {
        node $i [expr ($i -1)*$dl] 0 0;
}

fix $StartNode  1 1 1 1 1 1;


set ColSecTag 1
# define MATERIAL properties
set Es [expr 193050];               # Steel Young's Modulus N/mm2
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
set matID 1
uniaxialMaterial Elastic $matID $Es;
set J [expr  0.02473958*$mm4]
set GJ [expr $Gs*$J]
set z0 [expr 0.64625474*$mm];
set y0 [expr -0.68720012*$mm];

section FiberAsym $ColSecTag $y0 $z0 -GJ $GJ {
    fiber [expr -0.6872*$mm]  [expr 0.6463*$mm] [expr 0.0625*$mm2] $matID
    fiber [expr -0.5864*$mm]  [expr 0.4175*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr -0.4857*$mm]  [expr 0.1887*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr -0.3849*$mm] [expr -0.0401*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr -0.2841*$mm] [expr -0.2689*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr -0.1834*$mm] [expr -0.4977*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr -0.0826*$mm] [expr -0.7265*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.0182*$mm] [expr -0.9553*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.1189*$mm] [expr -1.1841*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.2197*$mm] [expr -1.4129*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.3205*$mm] [expr -1.6417*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.4212*$mm] [expr -1.8705*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  -0.458*$mm] [expr 0.7470*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  -0.229*$mm] [expr 0.8478*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  -0.0008*$mm] [expr 0.9486*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.2280*$mm]  [expr 1.0493*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.4568 *$mm] [expr 1.1501*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.6856*$mm]  [expr 1.2509*$mm]  [expr 0.0625*$mm2]  $matID
    fiber [expr  0.9143*$mm]  [expr 1.3516 *$mm] [expr 0.0625*$mm2]  $matID
}

set IDColTransf 1; # all members
set ColTransfType Corotational;           # options for columns: Linear PDelta Corotational
geomTransf $ColTransfType  $IDColTransf 0 0 1;    #define geometric transformation: performs a corotational geometric transformation
set numIntgrPts 2;        # number of Gauss integration points
    for {set i 1} {$i<$EndNode} {incr i 1} {
    set elemID $i
    set nodeI $i
    set nodeJ [expr $i+1]
    element mixedBeamColumnAsym $elemID $nodeI $nodeJ $numIntgrPts $ColSecTag $IDColTransf -shearCenter $y0 $z0;
}

# Define RECORDERS -------------------------------------------------------------
recorder Node -file DispMB6-$nl.out -time -node $EndNode -dof 1 2 3 4 5 6 disp;                      # displacements of middle node


# define Load-------------------------------------------------------------
set N 1000.0;
pattern Plain 2 Linear {
  # NodeID, Fx, Fy, Fz, Mx, My, Mz
  load $EndNode -$N 0 0 0 0 0;
}

# define ANALYSIS PARAMETERS------------------------------------------------------------------------------------
numberer Plain
system BandGeneral
constraints Plain
test NormUnbalance 1e-4 500 2
algorithm Newton;
analysis Static;

#define EQPath_Method_MRD 1 
#define EQPath_Method_NP 2
#define EQPath_Method_UNP 3
#define EQPath_Method_CYL 4
#define EQPath_Method_MNP 5
#define EQPath_Method_GDC 6
#define EQPath_Method_MUNP 7
#define EQPath_Method_GMRD 9
#define EQPath_Method_PEP 10

integrator EQPath 5 9 element -dof 4 5 6 

set i 0;
set ret 0;
set steps  10000;

while {$ret==0 && [getTime] < 200} {
	incr i
	set ret [analyze 1];
  #puts "$i [expr [getTime]] [expr [nodeDisp $EndNode 2]]"
}

puts "Finished"
#--------------------------------------------------------------------------------
set finishTime [clock clicks -milliseconds];
puts "Time taken: [expr ($finishTime-$startTime)/1000] sec"
set systemTime [clock seconds]
puts "Finished Analysis: [clock format $systemTime -format "%d-%b-%Y %H:%M:%S"]"