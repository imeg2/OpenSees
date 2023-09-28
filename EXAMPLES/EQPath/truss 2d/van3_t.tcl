
wipe

model basic -ndm 2 -ndf 2;

node 0 -9.6593 0  
node 1 0 2.5882  
node 2 9.6593 0  
node 3 0 0


fix 0 1 1
fix 2 1 1
fix 3 1 1


#element corotTruss $eleTag $iNode $jNode $secTag
uniaxialMaterial Elastic 1 2.0e5
element corotTruss 1 0 1 1 1 
element corotTruss 2 3 1 1 1
element corotTruss 3 2 1 1 1

pattern Plain 1 Linear { 
	load 1 1000 -1000; 
}

set targetNode 1
recorder Node -file van3.txt -time -node $targetNode -dof 1 2 disp

numberer Plain

system BandGeneral

constraints Transformation

test NormDispIncr 1E-07 25

algorithm Newton;

set arc_length 0.25;
set type 1;
integrator EQPath $arc_length $type

analysis Static;

for {set i 1} {$i < 2000} {incr i 1} {
	analyze 1
	if {[nodeDisp 1 2] < -6 || [nodeDisp 1 2] > 0.5} {
        puts "reach targe, exit"
      break
    }
}

system BandGeneral