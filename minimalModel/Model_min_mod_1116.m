function Func = Model_min_mod_1116(t,y,int,rdf,Dtot,Kir,Kr1,Kr2,Ksyn,Ksynr,KbI,Klri);
% solving ODEs
Func = zeros(3, 1);

% y(1)  BPtot             
% y(2)  LR-int4 first
% y(3)  BP-int4-rdf2 first

kpr=6; 
kmr1=kpr/Kr1;
kmr2=kpr/Kr2;
kpsyn=0.006;
kmsyn=kpsyn/Ksyn;
kpsynr=0.06; 
kmsynr=kpsynr/Ksynr;

intrdf=int*rdf/Kir;
Bp=(y(1)-y(3))/(1+int^4/KbI+intrdf^4/KbI+int^2*intrdf^2/KbI);
Lr=(Dtot-y(1)-y(2))/(1+int^4/KbI+intrdf^4/KbI+int^2*intrdf^2/Klri);
BpI=Bp*int^4/KbI;
LrI2=Lr*int^4/KbI;
LrIR=Lr*intrdf^4/KbI;
BpIR2=Bp*intrdf^4/KbI;

Func(1) = kmr1*y(2)-kpr*BpI+kpr*LrIR-kmr2*y(3);
Func(2) = kpr*BpI-kmr1*y(2)-kpsyn*y(2)+kmsyn*LrI2;
Func(3) = kpr*LrIR-kmr2*y(3)-kpsynr*y(3)+kmsynr*BpIR2;

