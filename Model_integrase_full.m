function Func = Model_integrase_full(t,y,Kii,Ki,Kir,Kbi1,Kbi2,Kbi3,Kbi4,Ks1,Ks2,Ks3,Ks4,Kr1,Kr2,Kmod,Kmodr,Ks01,Ks02);
 

Func = zeros(36, 1);

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

% kp - association rate; km - dissociation rate

kps01=40;
kps02=kps01;
kp=60;
kmii=Kii*kp;
kmir=Kir*kp;
kmbi1=Kbi1*kp;
kpbi2=10;
kmbi2=Kbi2*kpbi2;
kmbi3=Kbi3*kp;
kmbi4=Kbi4*kp;
kpi=3;
kmi=Ki*kpi;
kps1=0.8;
kps3=kps1;
kms01=Ks01*kps01;
kms02=Ks02*kps02;
kms1=Ks1*kps1;
kms3=Ks3*kps3;
kps2=0.0001;
kms2=Ks2*kps2;
kps4=0.005;
kms4=Ks4*kps4;
kpr=1; 
kmr1=kpr/Kr1;
kmr2=kpr/Kr2;
kpmod=1; 
kpmodr=kpmod; 
kmmod=kpmod/Kmod;
kmmodr=kpmodr/Kmodr;

Func(1) = kmbi2*y(6)-kpbi2*y(1)*y(4)+kmbi3*(y(13)+y(22))-kp*y(1)*(y(3)+y(20));
Func(2) = kmii*(2*y(4)+y(20))+kmir*y(19)-kp*(2*y(2)*y(2)+y(2)*y(19)+y(2)*y(23));
Func(3) = kp*(y(20)*y(23)+y(19)*y(19))-(kmir+kmii)*y(3)-kp*y(3)*(y(24)+y(12)+y(1)+y(13))-kps01*y(3)*(y(6)+y(22))-kps02*y(3)*(y(5)+y(21))+kmbi4*(y(12)+y(14))+kms01*(y(29)+y(30))+kms02*(y(26)+y(27))+kmbi3*(y(13)+y(15));
Func(4) = kp*y(2)*y(2)-kmii*y(4)-kp*y(4)*(y(23)+y(24)+y(5))-kpbi2*y(4)*(y(1)+y(6))-kps01*y(4)*(y(13)+y(22))-kps02*y(4)*(y(12)+y(21))+kmir*y(20)+kmbi1*(y(5)+y(7))+kms01*(y(28)+y(29))+kms02*(y(25)+y(26))+kmbi2*(y(6)+y(8))-(kpi*y(4)*(y(7)+y(14)+y(25)+y(26)+y(27))-kmi*(y(31)+y(32)+y(34)+y(35)+y(36)));
Func(5) = kp*y(24)*y(4)-kmbi1*y(5)-kp*y(5)*y(4)-kps02*y(5)*(y(20)+y(3))+kmbi1*y(7)+kms02*(y(25)+y(26))-kp*y(23)*y(5)+kmir*y(21);
Func(6) = kpbi2*y(4)*(y(1)-y(6))-kmbi2*(y(6)-y(8))-kps01*y(6)*(y(20)+y(3))+kms01*(y(28)+y(29))-kp*y(23)*y(6)+kmir*y(22);
Func(7) = kp*y(5)*y(4)-kmbi1*y(7)-kps1*y(7)+kms1*y(10)-(kpi*y(4)*y(7)-y(31)*kmi);
Func(8) = kpbi2*y(6)*y(4)-kmbi2*y(8)-kps2*y(8)+kms2*y(9);
Func(9) = kpmod*y(11)-kmmod*y(9)-kms2*y(9)+kps2*y(8);
Func(10) = kps1*y(7)-kms1*y(10)-kpr*y(10)+kmr1*y(11);
Func(11) = kpr*y(10)-kmr1*y(11)-kpmod*y(11)+kmmod*y(9);
Func(12) = kp*y(24)*y(3)-kmbi4*y(12)+kp*y(21)*y(23)-kmir*y(12)-kp*y(12)*y(3)-kps02*y(12)*(y(4)+y(20))+kms02*(y(26)+y(27))+kmbi4*y(14);
Func(13) = kp*y(1)*y(3)-kmbi3*y(13)+kp*y(22)*y(23)-kmir*y(13)-kp*y(13)*y(3)-kps01*y(13)*(y(4)+y(20))+kms01*(y(29)+y(30))+kmbi3*y(15);
Func(14) = kp*y(3)*y(12)-kmbi4*y(14)+kms4*y(16)-kps4*y(14)-(kpi*y(4)*y(14)-y(32)*kmi);
Func(15) = kp*y(3)*y(13)-kmbi3*y(15)-kps3*y(15)+kms3*y(18);
Func(16) = kpmodr*y(17)-kmmodr*y(16)-kms4*y(16)+kps4*y(14);
Func(17) = kpr*y(18)-kmr2*y(17)+kmmodr*y(16)-kpmodr*y(17);
Func(18) = kps3*y(15)-kms3*y(18)-kpr*y(18)+kmr2*y(17);
Func(19) = kp*(y(2)*y(23)-y(2)*y(19)-2*y(19)*y(19))-kmir*y(19)+kmii*(y(20)+2*y(3));
Func(20) = kp*(y(4)*y(23)+y(2)*y(19))-(kmir+kmii)*y(20)-kp*y(23)*y(20)+kmir*y(3)-kp*y(20)*(y(24)+y(1))-kps01*y(20)*(y(6)+y(22)+y(13))-kps02*y(20)*(y(5)+y(21)+y(12))+kmbi4*y(21)+kms01*(y(28)+y(29)+y(30))+kms02*(y(25)+y(26)+y(27))+kmbi3*y(22);
Func(21) = kp*y(24)*y(20)-kmbi4*y(21)-kps02*y(21)*(y(4)+y(20)+y(3))+kms02*(y(25)+y(26)+y(27))+kp*y(23)*y(5)-kmir*y(21)-kp*y(21)*y(23)+kmir*y(12);
Func(22) = kp*y(1)*y(20)-kmbi3*y(22)-kps01*y(22)*(y(20)+y(3)+y(4))+kms01*(y(28)+y(29)+y(30))+kp*y(23)*y(6)-kmir*y(22)-kp*y(22)*y(23)+kmir*y(13);
Func(23) = kmir*(y(19)+y(20)+y(3)+y(21)+y(22)+y(12)+y(13))-kp*y(23)*(y(2)+y(4)+y(20)+y(5)+y(6)+y(21)+y(22));
Func(24) = kmbi1*y(5)+kmbi4*(y(12)+y(21))-kp*y(24)*(y(4)+y(20)+y(3));
Func(25) = kps02*(y(4)*y(21)+y(20)*y(5))-y(25)*2*kms02-(kpi*y(4)*y(25)-y(34)*kmi);
Func(26) = kps02*(y(4)*y(12)+y(20)*y(21)+y(3)*y(5))-y(26)*3*kms02-(kpi*y(4)*y(26)-y(35)*kmi);
Func(27) = kps02*(y(20)*y(12)+y(3)*y(21))-y(27)*2*kms02-(kpi*y(4)*y(27)-y(36)*kmi);
Func(28) = kps01*(y(4)*y(22)+y(20)*y(6))-y(28)*2*kms01;
Func(29) = kps01*(y(4)*y(13)+y(20)*y(22)+y(3)*y(6))-y(29)*3*kms01;
Func(30) = kps01*(y(20)*y(13)+y(3)*y(22))-y(30)*2*kms01;
Func(31) = kpi*y(4)*y(7)-y(31)*kmi;
Func(32) = kpi*y(4)*y(14)-y(32)*kmi;
Func(33) = 0;
Func(34) = kpi*y(4)*y(25)-y(34)*kmi;
Func(35) = kpi*y(4)*y(26)-y(35)*kmi;
Func(36) = kpi*y(4)*y(27)-y(36)*kmi;

