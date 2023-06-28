load("results_Tachidius_discipes.mat")

% compute temperature correction factors
pars_T = T_A;
TC_ah = tempcorr(temp.ah, T_ref, pars_T);
TC_tL = tempcorr(temp.tL, T_ref, pars_T);
TC_tL_18 = tempcorr(temp.tL_18, T_ref, pars_T);

metaPar.model = 'abp'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.T_A = 241.503*TC_tL;