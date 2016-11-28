% run this file to reproduce Fig.4
% effect of addition of RDF on PxB(-R) reaction. RDF added after 1 hour (solid line)
% The kinetics of the LxR(+R) reaction is shown for comparison by a dotted line
% The computations were performed with 400 nM integrase and 800 nM RDF

Dtot=0.01; Kr1=2.8; Kr2=2; Kir=0.05; Ksyn=0.36; Ksynr=0.5; KbI=0.0001; Klri=0.00002; int_tot=0.4;

y0=[Dtot 0 0]; 
rdf_tot=0; 
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1);
PBt=Y(:,1); 

figure()
p1=semilogx(T,PBt/Dtot,'b','DisplayName','PxB(-R) reaction');
hold on;

rdf_tot=0.8;
y0=[Y(end,1) Y(end,2) Y(end,3)]; 
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 100);
PBt=Y(:,1); 
T=T+1;

semilogx(T,PBt/Dtot,'b');

y0=[0 0 0]; 
rdf_tot=0; 
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1);
PBt=Y(:,1); 

p2=semilogx(T,PBt/Dtot,'b:','DisplayName','LxR(+R) reaction');

rdf_tot=0.8;
y0=[Y(end,1) Y(end,2) Y(end,3)]; 
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 100);
PBt=Y(:,1); 
T=T+1;

semilogx(T,PBt/Dtot,'b:');
xlabel('time, h');
ylabel('PB, fraction from total DNA');
legend([p1,p2],'PxB(-R) reaction','LxR(+R) reaction')
title('Simulated effect of addition of RDF after 1h')
