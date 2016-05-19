Kr1=1; Kr2=1; Dtot=0.01; Ki=0.02; Kii=0.3; Kir=0.05; Kmod=3.4; Kmodr=1.9; Ks01=0.001; Ks02=0.007;
Ks1=0.1; Ks2=0.12; Ks3=0.1; Ks4=0.013; Kbi1=0.02; Kbi2=0.01; Kbi3=0.025; Kbi4=0.05; %Ks1,Ks3,Kb1,Kb3-dissociation constants

% time unit - min

% y(1)   LR                   
% y(2)   int
% y(3)   int2-rdf2
% y(4)   int2
% y(5)   BP-int2
% y(6)   LR-int2
% y(7)   BP-int4
% y(8)   LR-int4
% y(9)   LR-int4 synapse second
% y(10)  BP-int4 synapse 
% y(11)  LR-int4 synapse first
% y(12)  BP-int2-rdf2
% y(13)  LR-int2-rdf2
% y(14)  BP-int4-rdf4
% y(15)  LR-int4-rdf4
% y(16)  BP-int4-rdf4 synapse second
% y(17)  BP-int4-rdf4 synapse first
% y(18)  LR-int4-rdf4 synapse
% y(19)  int-rdf
% y(20)  int2-rdf
% y(21)  BP-int2-rdf
% y(22)  LR-int2-rdf
% y(23)  rdf
% y(24)  BP
% y(25)  BP-int4-rdf
% y(26)  BP-int4-rdf2
% y(27)  BP-int4-rdf3
% y(28)  LR-int4-rdf
% y(29)  LR-int4-rdf2
% y(30)  LR-int4-rdf3
% y(31)  BP-int6i
% y(32)  BP-int6-rdf4i
% y(33)  
% y(34)  BP-int6-rdfi
% y(35)  BP-int6-rdf2i
% y(36)  BP-int6-rdf3i

y0=[0 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0. Dtot 0 0 0 0 0 0 0 0 0 0 0 0]; % initial conditions for PxB(-RDF)
%y0=[Dtot 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.8 0 0 0 0 0 0 0 0 0 0 0 0 0]; % initial conditions for LxR(-RDF)

tJ=[0 1 2 4 8 16 32 64 120 180]; % kinetics data
prJ1=[0 28 36 47 56 64 68 72 73 74]; % 0.4 int + BP
pr_err1=[0 2 1 1 2 2 1 2 2 2];
prJ2=[0 24 32 42 51 57 61 64 66 67]; % 0.4 int + 0.8 rdf + LR
pr_err2=[0 3 2 1 1 1 1 1 2 2];

options = odeset(); 

    t=[0 180];
    [T, Y] = ode15s(@Model_integrase_full,t,y0,options,Kii,Ki,Kir,Kbi1,Kbi2,Kbi3,Kbi4,Ks1,Ks2,Ks3,Ks4,Kr1,Kr2,Kmod,Kmodr,Ks01,Ks02);

LRt=Y(:,1)+Y(:,6)+Y(:,8)+Y(:,9)+Y(:,11)+Y(:,13)+Y(:,15)+Y(:,18)+Y(:,22)+Y(:,28)+Y(:,29)+Y(:,30);
BPt=Y(:,5)+Y(:,7)+Y(:,10)+Y(:,12)+Y(:,14)+Y(:,16)+Y(:,17)+Y(:,21)+Y(:,24)+Y(:,25)+Y(:,26)+Y(:,27)+Y(:,31)+Y(:,32)+Y(:,34)+Y(:,35)+Y(:,36);
int_tot=Y(:,2)+Y(:,19)+2*(Y(:,3)+Y(:,20)+Y(:,21)+Y(:,22)+Y(:,4)+Y(:,5)+Y(:,6)+Y(:,12)+Y(:,13))+4*(Y(:,7)+Y(:,8)+Y(:,9)+Y(:,10)+Y(:,11)+Y(:,14)+Y(:,15)+Y(:,16)+Y(:,17)+Y(:,18)+Y(:,28)+Y(:,29)+Y(:,30)+Y(:,25)+Y(:,26)+Y(:,27))+6*(Y(:,31)+Y(:,32)+Y(:,34)+Y(:,35)+Y(:,36));
rdf_tot=Y(:,23)+Y(:,19)+Y(:,20)+Y(:,21)+Y(:,22)+Y(:,25)+Y(:,28)+Y(:,34)+2*(Y(:,3)+Y(:,12)+Y(:,13)+Y(:,26)+Y(:,29)+Y(:,35))+3*(Y(:,27)+Y(:,30)+Y(:,36))+4*(Y(:,14)+Y(:,15)+Y(:,16)+Y(:,17)+Y(:,18)+Y(:,32));
DNA_rdf=Y(:,13)+Y(:,15)+Y(:,18)+Y(:,22)+Y(:,28)+Y(:,29)+Y(:,30)+Y(:,12)+Y(:,14)+Y(:,16)+Y(:,17)+Y(:,21)+Y(:,25)+Y(:,26)+Y(:,27)+Y(:,32)+Y(:,34)+Y(:,36);

figure (1)
plot(T,LRt/Dtot,'r');
hold on;
plot(T,BPt/Dtot,'b');
hold on;
plot(T,Y(:,9)/Dtot,'r--');
hold on;
plot(T,Y(:,11)/Dtot,'r:');
hold on;
plot(T,Y(:,16)/Dtot,'b--');
hold on;
plot(T,Y(:,17)/Dtot,'b:');
hold on;
title({'LRtot-red (LR-int_s_2 - dash; LR-int_s_1 - dot)';'BPtot-blue (BP-int-rdf_s_2- dash; BP-int-rdf_s_1- dot)'});

figure (2)
plot(T,BPt/Dtot*100,'b');
hold on;
plot(T,LRt/Dtot*100,'r');
hold on;
errorbar(tJ,prJ1,pr_err1,'r:');
hold on;
errorbar(tJ,prJ2,pr_err2,'b:');
hold on;
title('LRtot(%)-red; BPtot(%)-blue');

