% run this file to reproduce Fig.2D
% Dependence of the level of PB product (% from max) after 3h of LxR reaction on concentration of integrase 
% Different lines correspond to different concentrations of RDF: 0 (blue),
% 50 nM (red), 100 nM (yellow), 200 nM (black), 400 nM (magenda), 800 nM (green)

Dtot=0.01; 

y0=[0 0 0]; % initial conditions for LxR reaction

RDF = [0,0.05,0.1,0.2,0.4,0.8];
col = ['b','r','y','k','m','g'];

figure()
hold on;
for r=1:length(RDF)


    rdf_tot=RDF(r); %concentrations of RDF
    for i=1:17
        x_int(i)=0.05*(i-1);
        int_tot=x_int(i);
        [Y,T]=min_mod_251116(rdf_tot, int_tot, y0, 3);
        PB=Y(:,1); 
        PB_3h(i)=PB(end)/Dtot*100;
    end

    plot(x_int*1000,PB_3h,col(r),'DisplayName',num2str(rdf_tot*1000));

end

xlabel('integrase, nM');
ylabel('PB, %');
title('PB product after 3h for different integrase and RDF concentrations (nM)')
legend('show')

