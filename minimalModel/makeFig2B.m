% run this file to reproduce Fig.2B
% Dependence of the level of LR product (% from max) after 3h of PxB reaction on concentration of integrase 
% Different lines correspond to different concentrations of RDF: 0 (blue),
% 50 nM (red), 100 nM (yellow), 200 nM (black), 400 nM (magenda), 800 nM (green)

Dtot=0.01; 

y0=[Dtot 0 0]; % initial conditions for PxB reaction

RDF = [0,0.05,0.1,0.2,0.4,0.8];
col = ['b','r','y','k','m','g'];

figure()
for r=1:length(RDF)
    
    rdf_tot=RDF(r); %concentrations of RDF
    %y0=[0 0 0]; rdf_tot=0.4; % initial conditions for LxR reaction
    for i=1:17
        x_int(i)=0.05*(i-1);
        int_tot=x_int(i);
        [Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 3);
        PB=Y(:,1); 
        LRt=Dtot-PB; 
        LR_3h(i)=LRt(end)/Dtot*100;
    end

    plot(x_int*1000,LR_3h,col(r),'DisplayName',num2str(rdf_tot*1000));
    hold on;

end

xlabel('integrase, nM');
ylabel('LR, %');
title('LR product after 3h for different integrase and RDF concentrations (nM)')
legend('show')
