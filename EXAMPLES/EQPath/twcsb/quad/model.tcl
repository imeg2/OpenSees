wipe;
model basic -ndm 3 -ndf 6;

source LibUnits.tcl;

source "nodes.tcl";

fix 1  1 1 1 1 1 1;

set StartNode 1;
set EndNode 7;
set ColSecTag 1
# define MATERIAL properties
set Es [expr 27910.0*$ksi];               # Steel Young's Modulus
set nu 0.3;
set Gs [expr $Es/2./[expr 1+$nu]];  # Torsional stiffness Modulus
set matID 1
uniaxialMaterial Elastic $matID $Es;
set J [expr  0.02473958*$in4]
set GJ [expr $Gs*$J]
set z0 [expr 0.64625474*$in];
set y0 [expr -0.68720012*$in];

section FiberAsym $ColSecTag $y0 $z0 -GJ $GJ {
	
}
section('FiberAsym', ColSecTag, y0, z0, '-GJ', Gj);
fiber -0.6872,  0.6463, 0.0625, matID);
fiber(-0.5864,  0.4175, 0.0625, matID);
fiber(-0.4857,  0.1887, 0.0625, matID);
fiber(-0.3849, -0.0401, 0.0625, matID);
fiber(-0.2841, -0.2689, 0.0625, matID);
fiber(-0.1834, -0.4977, 0.0625, matID);
fiber(-0.0826, -0.7265, 0.0625, matID);
fiber( 0.0182, -0.9553, 0.0625, matID);
fiber( 0.1189, -1.1841, 0.0625, matID);
fiber( 0.2197, -1.4129, 0.0625, matID);
fiber( 0.3205, -1.6417, 0.0625, matID);
fiber( 0.4212, -1.8705, 0.0625, matID);
fiber( -0.4584, 0.7470, 0.0625, matID);
fiber( -0.2296, 0.8478, 0.0625, matID);
fiber( -0.0008, 0.9486, 0.0625, matID);
fiber( 0.2280,  1.0493, 0.0625, matID);
fiber( 0.4568,  1.1501, 0.0625, matID);
fiber( 0.6856,  1.2509, 0.0625, matID);
fiber( 0.9143,  1.3516, 0.0625, matID);


set targetNode 543
#set targetNode 1805

pattern Plain 1 Linear { 
    load $targetNode 0 -1 0 0 0 0
}

numberer Plain

system SparseGeneral 

constraints Plain

test NormUnbalance 1e-5 20 2

algorithm KrylovNewton;

#define EQPath_Method_MRD 1 
#define EQPath_Method_NP 2
#define EQPath_Method_UNP 3
#define EQPath_Method_CYL 4
#define EQPath_Method_MNP 5
#define EQPath_Method_GDC 6
#define EQPath_Method_MUNP 7
#define EQPath_Method_GMRD 9
#define EQPath_Method_PEP 10

set arc_length 5;
set type 9;

recorder Node -file node-$targetNode-method-$type-KrylovNewton-element.txt  -time -node $targetNode -dof 1 2 3 disp

integrator EQPath $arc_length $type element -dof 4 5 6

analysis Static;

puts "start"

set i 0;
set ret 0;
set steps  1600;
while {$i<$steps && $ret==0} {
	incr i
	set ret [analyze 1];
	if {[getTime] > 20 } {
      puts "reach targe, exit"
      break
    }

puts "[getTime]"
}

puts "done"
remove recorders