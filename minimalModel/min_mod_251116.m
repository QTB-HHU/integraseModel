function [Y,T] = min_mod_251116(rdf_tot, int_tot, y0,Tfin)
Dtot=0.01; Kr1=2.8; Kr2=2; Kir=0.05; Ksyn=0.36; Ksynr=0.5; KbI=0.0001; Klri=0.00002;

% time unit - hour

%y0 - initial conditions for the reaction
%int_tot, rdf_tot - concentrations of integrase and RDF in mkM

options = odeset(); 
    %options = odeset('MaxStep',0.0001);
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;

    t=[0 Tfin]; % time interval
    % kinetics of the reaction
        options = odeset('MaxStep',0.1);
    [T, Y] = ode15s(@Model_min_mod_1116,t,y0,options,int,rdf,Dtot,Kir,Kr1,Kr2,Ksyn,Ksynr,KbI,Klri);

PB=Y(:,1); % PB total


 

