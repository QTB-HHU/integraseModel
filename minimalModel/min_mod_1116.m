% main model file, calling for ode file Model_min_mod_1116.m

Kr1=2.8; Kr2=2; Kir=0.05; Ksyn=0.36; Ksynr=0.5; KbI=0.0001; Klri=0.00002;

% time unit - hour

% y(1)  BPtot             
% y(2)  LR-int4 first
% y(3)  BP-int4-rdf2 first

Dtot=0.01; y0=[Dtot 0 0]; % initial conditions for PxB reaction
int_tot=0.4; rdf_tot=0.; %concentrations of integrase and RDF in mkM
%y0=[0 0 0]; rdf_tot=0.4; % initial conditions for LxR reaction

options = odeset(); 
    %options = odeset('MaxStep',0.0001);
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;

    t=[0 3]; % time interval
    [T, Y] = ode15s(@Model_min_mod_1116,t,y0,options,int,rdf,Dtot,Kir,Kr1,Kr2,Ksyn,Ksynr,KbI,Klri);

LRt=Dtot-Y(:,1);
Intrdf=int*rdf/Kir;
BP=(Y(:,1)-Y(:,3))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/KbI);
LR=(Dtot-Y(:,1)-Y(:,2))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/Klri);
BPI=BP*int^4/KbI;
LRI2=LR*int^4/KbI;
LRIR=LR*Intrdf^4/KbI;
BPIR2=BP*Intrdf^4/KbI;

% kinetics of the total LR and PB during PxB reaction; t=[0 3]; integrase=0.4 mkM, rdf=0
% to run LxR reaction, change initial condition to y0=[0 0 0]
% to calculate the product level at 3h, use the last datapoint of the product vector: LRt for PxB reaction; Y(:1) for LxR reaction
% to calculate reaction kinetics with other integrase or RDF concentration, change values of variables int_tot=0.4 and rdf_tot

figure (1)
plot(T,LRt/Dtot,'r');
hold on;
plot(T,Y(:,1)/Dtot,'b');
hold on;
title('LR_t-red; BP_t-blue');
 

