% run this file to reproduce Fig.3
% time courses of the most abundant DNA-containing products
% Panels A and B show the kinetics of the “allowed” reactions (PxB(-R) and LxR(+R))
% Panels C and D display the kinetics of the “forbidden” reactions (LxR(-R) and PxB(+R))

Dtot=0.01; Kr1=2.8; Kr2=2; Kir=0.05; Ksyn=0.36; Ksynr=0.5; KbI=0.0001; Klri=0.00002; int_tot=0.4;

figure()
% panel A
subplot(2,2,1)

y0=[Dtot 0 0]; % initial conditions for PxB reaction
rdf_tot=0; %concentrations of RDF
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;
Intrdf=int*rdf/Kir;
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1000);
PB=Y(:,1); 
LRt=(Dtot-PB); 
LR=(Dtot-Y(:,1)-Y(:,2))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/Klri);
LRI2=LR*int^4/KbI;

semilogx(T,LRt/Dtot,'r');
hold on;
semilogx(T,Y(:,2)/Dtot,'r--');
semilogx(T,LRI2/Dtot,'r:');
xlim([1e-3,1e3])
ylim([0,1])
ax=gca
ax.set('XTick',10.^[-3:3])
xlabel('time, h');
ylabel('LR, fraction from total DNA');
legend('LR_{tot}','LRI_1','LRI_2','Location','northwest')
title('"allowed" reaction PxB(-R)')

% panel B
y0=[0 0 0]; % initial conditions for LxR reaction
rdf_tot=0.8; %concentrations of RDF
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;
Intrdf=int*rdf/Kir;
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1000);
PBt=Y(:,1); 
BP=(Y(:,1)-Y(:,3))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/KbI);
BPIR2=BP*Intrdf^4/KbI;

subplot(2,2,2)
semilogx(T,PBt/Dtot,'b');
hold on;
semilogx(T,Y(:,3)/Dtot,'b--');
semilogx(T,BPIR2/Dtot,'b:');
xlim([1e-3,1e3])
ylim([0,1])
ax=gca;
ax.set('XTick',10.^[-3:3])
xlabel('time, h');
ylabel('PB, fraction from total DNA');
legend('PB_{tot}','PBIR_1','PBIR_2','Location','northwest')
title('"allowed" reaction LxR(+R)')

% panel C
y0=[0 0 0]; % initial conditions for LxR reaction
rdf_tot=0.; %concentrations of RDF
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;
Intrdf=int*rdf/Kir;
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1000);
PB=Y(:,1); 
LRt=(Dtot-PB); 
LR=(Dtot-Y(:,1)-Y(:,2))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/Klri);
LRI2=LR*int^4/KbI;

subplot(2,2,3)
semilogx(T,LRt/Dtot,'r');
hold on;
semilogx(T,Y(:,2)/Dtot,'r--');
semilogx(T,LRI2/Dtot,'r:');
xlim([1e-3,1e3])
ylim([0,1])
ax=gca;
ax.set('XTick',10.^[-3:3])
xlabel('time, h');
ylabel('LR, fraction from total DNA');
legend('LR_{tot}','LRI_1','LRI_2','Location','northwest')
title('"forbidden" reaction LxR(-R)')

% panel D
y0=[Dtot 0 0]; % initial conditions for PxB reaction
rdf_tot=0.8; %concentrations of RDF
b=rdf_tot-int_tot+Kir;
int=0.5*(sqrt(b*b+4*int_tot*Kir)-b);
rdf=rdf_tot-int_tot+int;
Intrdf=int*rdf/Kir;
[Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 1000);
PBt=Y(:,1); 
BP=(Y(:,1)-Y(:,3))/(1+int^4/KbI+Intrdf^4/KbI+int^2*Intrdf^2/KbI);
BPIR2=BP*Intrdf^4/KbI;

subplot(2,2,4)
semilogx(T,PBt/Dtot,'b');
hold on;
semilogx(T,Y(:,3)/Dtot,'b--');
semilogx(T,BPIR2/Dtot,'b:');
xlim([1e-3,1e3])
ylim([0,1])
ax=gca;
ax.set('XTick',10.^[-3:3])
xlabel('time, h');
ylabel('PB, fraction from total DNA');
legend('PB_{tot}','PBIR_1','PBIR_2','Location','northwest')
title('"forbidden" reaction PxB(+R)')

