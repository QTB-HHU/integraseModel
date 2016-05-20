# integraseModel
MATLAB files to simulate in vitro experiments of Integrase with or without Reaction Directionality Factors (RDFs)

To simulate reaction kinetics for 3 hours (Figs. 5-9 in the manuscript), execute the file "integrase_full.m" in MATLAB.

To simulate the different reactions, use the following initial conditions (for 0.4 µM integrase, 0.8 µM RDF) in line 43 of the file "integrase_full.m":

1) BxP(-RDF)
y0=[0 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0. Dtot 0 0 0 0 0 0 0 0 0 0 0 0]; 

2) LxR(+RDF)
y0=[Dtot 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.8 0 0 0 0 0 0 0 0 0 0 0 0 0];

3) BxP(+RDF)
y0=[0 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0.8 Dtot 0 0 0 0 0 0 0 0 0 0 0 0];

4) LxR(-RDF)
y0=[Dtot 0.4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0. 0 0 0 0 0 0 0 0 0 0 0 0 0];

Execution of "integrase_full.m" produces two plots:
Plot 2 will show the kinetics of the total amounts of LR and BP (% of the total DNA concentration Dtot), together with experimental data. Plot 1 will show the kinetics of major fractions of the DNA product, together with the total amounts of LR and BP, normalized to Dtot.

To simulate Fig. 9A use initial condition 1 and increase 5-fold the value of Kr1 or Kmod and in the same time decrease 5-fold the value of Ks2 to keep energy conserved. To simulate Fig. 9B use initial condition 2 and increase 5-fold the value of Kr2 or Kmodr and in the same time decrease 5-fold the value of Ks4.
