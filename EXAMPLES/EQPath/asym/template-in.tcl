set systemTime [clock seconds]
puts "Starting Analysis: [clock format $systemTime -format "%d-%b-%Y %H:%M:%S"]"
set startTime [clock clicks -milliseconds];
# SET UP ----------------------------------------------------------------------------
wipe;                             # clear memory of all past model definitions
model BasicBuilder -ndm 3 -ndf 6; # Define the model builder, ndm=#dimension, ndf=#dofs

source LibUnits.tcl;                      # define units
set length [expr 1400/$mm]
set nl {nl}
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
set Es [expr 27910.0];               # Steel Young's Modulus N/mm2
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
set matID 1
uniaxialMaterial Elastic $matID $Es;
set J [expr  0.02473958]
set GJ [expr $Gs*$J]
set z0 [expr 0.64625474];
set y0 [expr -0.68720012];

section FiberAsym $ColSecTag $y0 $z0 -GJ $GJ {
    fiber [expr -0.6872*1]  [expr 0.6463*1] [expr 0.0625*1] $matID
    fiber [expr -0.5864*1]  [expr 0.4175*1]  [expr 0.0625*1]  $matID
    fiber [expr -0.4857*1]  [expr 0.1887*1]  [expr 0.0625*1]  $matID
    fiber [expr -0.3849*1] [expr -0.0401*1]  [expr 0.0625*1]  $matID
    fiber [expr -0.2841*1] [expr -0.2689*1]  [expr 0.0625*1]  $matID
    fiber [expr -0.1834*1] [expr -0.4977*1]  [expr 0.0625*1]  $matID
    fiber [expr -0.0826*1] [expr -0.7265*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.0182*1] [expr -0.9553*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.1189*1] [expr -1.1841*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.2197*1] [expr -1.4129*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.3205*1] [expr -1.6417*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.4212*1] [expr -1.8705*1]  [expr 0.0625*1]  $matID
    fiber [expr  -0.458*1] [expr 0.7470*1]  [expr 0.0625*1]  $matID
    fiber [expr  -0.229*1] [expr 0.8478*1]  [expr 0.0625*1]  $matID
    fiber [expr  -0.0008*1] [expr 0.9486*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.2280*1]  [expr 1.0493*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.4568 *1] [expr 1.1501*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.6856*1]  [expr 1.2509*1]  [expr 0.0625*1]  $matID
    fiber [expr  0.9143*1]  [expr 1.3516 *1] [expr 0.0625*1]  $matID
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


# define Load-------------------------------------------------------------
set N 224.809;
pattern Plain 2 Linear {
  # NodeID, Fx, Fy, Fz, Mx, My, Mz
  load $EndNode -$N 0 0 0 0 0;
}

# define ANALYSIS PARAMETERS------------------------------------------------------------------------------------
numberer Plain
system BandGeneral
constraints Plain
test NormUnbalance 1e-4 100 2
algorithm Newton  ;
analysis Static   ;# define type of analysis static or transient

recorder Node -file [pwd]\\outputs\\{nl}\\{filename}.txt  -time -node $EndNode -dof 1 disp
integrator EQPath {arc_length} {type} {over} {dofs}

set i 0;
set ret 0;
set steps  10000;
set fileID [open [pwd]\\outputs\\{nl}\\{filename}.iter w]

while {$ret==0 && [getTime] < 100} {
	incr i
	set ret [analyze 1];
    puts  $fileID   "$i	[getTime]	[nodeDisp $EndNode 1]	[getCurrentIterations]	$ret"
	puts "[expr [getTime]]"
}

flush  $fileID       // flush any output that has been buffered    
close  $fileID       // close the file