wipe
model basic -ndm 2 -ndf 3;

node 0 50 0
node 1 49.6354437 6.026834013
node 2 48.54709087 11.96578321
node 3 46.75081213 17.73024435
node 4 44.27280128 23.2361586
node 5 41.14919329 28.40323734
node 6 37.42553741 33.15613291
node 7 33.15613291 37.42553741
node 8 28.40323734 41.14919329
node 9 23.2361586 44.27280128
node 10 17.73024435 46.75081213
node 11 11.96578321 48.54709087
node 12 6.026834013 49.6354437
node 13 3.06287E-15 50
node 14 -6.026834013 49.6354437
node 15 -11.96578321 48.54709087
node 16 -17.73024435 46.75081213
node 17 -23.2361586 44.27280128
node 18 -28.40323734 41.14919329
node 19 -33.15613291 37.42553741
node 20 -37.42553741 33.15613291
node 21 -41.14919329 28.40323734
node 22 -44.27280128 23.2361586
node 23 -46.75081213 17.73024435
node 24 -48.54709087 11.96578321
node 25 -49.6354437 6.026834013
node 26 -50 6.12574E-15


fix 0 1 1 0
fix 26 1 1 0


element	elastic2dGNL	0	0	1	10	200	1
element	elastic2dGNL	1	1	2	10	200	1
element	elastic2dGNL	2	2	3	10	200	1
element	elastic2dGNL	3	3	4	10	200	1
element	elastic2dGNL	4	4	5	10	200	1
element	elastic2dGNL	5	5	6	10	200	1
element	elastic2dGNL	6	6	7	10	200	1
element	elastic2dGNL	7	7	8	10	200	1
element	elastic2dGNL	8	8	9	10	200	1
element	elastic2dGNL	9	9	10	10	200	1
element	elastic2dGNL	10	10	11	10	200	1
element	elastic2dGNL	11	11	12	10	200	1
element	elastic2dGNL	12	12	13	10	200	1
element	elastic2dGNL	13	13	14	10	200	1
element	elastic2dGNL	14	14	15	10	200	1
element	elastic2dGNL	15	15	16	10	200	1
element	elastic2dGNL	16	16	17	10	200	1
element	elastic2dGNL	17	17	18	10	200	1
element	elastic2dGNL	18	18	19	10	200	1
element	elastic2dGNL	19	19	20	10	200	1
element	elastic2dGNL	20	20	21	10	200	1
element	elastic2dGNL	21	21	22	10	200	1
element	elastic2dGNL	22	22	23	10	200	1
element	elastic2dGNL	23	23	24	10	200	1
element	elastic2dGNL	24	24	25	10	200	1
element	elastic2dGNL	25	25	26	10	200	1





set targetNode 14
pattern Plain 1 Linear { 
    load $targetNode 0 -1 0
}

numberer Plain
system SparseSPD
constraints Transformation
algorithm Newton;

#define EQPath_Method_MRD 1
#define EQPath_Method_NP 2
#define EQPath_Method_UNP 3
#define EQPath_Method_CYL 4
#define EQPath_Method_MNP 5
#define EQPath_Method_GDC 6
#define EQPath_Method_MUNP 7
#define EQPath_Method_GMRD 9
#define EQPath_Method_PEP 10


set arc_length 15;
set type 9;
set error 1e-5;

recorder Node -file node-14-method-$type.rec  -time -node $targetNode -dof 2 disp

test NormDispIncr $error 20 2

integrator EQPath $arc_length $type element -dof 1 2 3

analysis Static;

puts "start"

set i 0;
set ret 0;
set steps  5000;
while {$i<$steps && $ret==0} {
	incr i
	set ret [analyze 1];
    if {[getTime] > 5.5 || [nodeDisp $targetNode 2] > 0} {
        puts "reach targe [getTime], exit"
      break
    }
}

puts "done"
remove recorders