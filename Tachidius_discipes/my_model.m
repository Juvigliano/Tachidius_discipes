%alalallalaa
% intialization has components
%1: reading the parameters from the matlab file. Depneding of what I prefer
%I can read it from the pars_init or.mat
load('results_Tachidius_discipes.mat','par') %this is to load the parameters from .mat
%into a structure called pars
% unpack par, data, auxData
cPar = parscomp_st(par);  vars_pull(par);
%2: a place where you specify food and temp.
T=297.15; %temperature = 24C
f=1;
%3: define the time you want to model
t_max=30;
t= linspace(0,t_max,1000);
%second part: think about what i want to predict
%i.e: length and weight. then I build a code for predicting length and
%weight for the time I specified (I can use code from my predict file).

[tsort, ~, ci] = unique(tL24(:,1)); 
[~, ELHR] = ode45(@dget_ELHR_sbp, tsort, ELHR0,[], p, TC24, f);
EL24 = ELHR(ci,2)/ del_M; 
prdData.tL24 = EL24;
% pack to output
%third: print out results. and save results into files

function dELHR = dget_ELHR_sbp(t, ELHR, p, TC, f)
  % Define changes in the state variables for abj model
  % t: time
  % ELHR: 4-vector with state variables
  %         E , J, reserve energy
  %         L , cm, structural length
  %         E_H , J , cumulated energy inversted into maturity (E_H in Kooijman 2010)
  %         E_R , J, reproduction buffer (E_R in Kooijman 2010)
  %         
  % dELHR: 4-vector with change in E, L, H, R

  % unpack state variables

  E = ELHR(1); L = ELHR(2); E_H = ELHR(3);

  % unpack par
  p_Am = p(1); v = p(2); p_M = p(3); k_J = p(4); 
  kap = p(5); kap_G = p(6); 
  E_G = p(7); E_Hb = p(8); E_Hp = p(9);

  % temp correction
  pT_Am = TC * p_Am ;
  vT = TC * v;  
  pT_M = TC * p_M;
  kT_J = TC * k_J; 

  pA = (pT_Am * f * L^2) * (E_H >= E_Hb);
  if E_H < E_Hp
    if  kap * E * vT >= pT_M * L^4 % section 4.1.5 comments to Kooy2010
        r = (E * vT/ L - pT_M * L^3/ kap)/ (E + E_G * L^3/ kap); % d^-1, specific growth rate  
    else 
        r = (E * vT/ L - pT_M * L^3/ kap)/ (E + kap_G * E_G * L^3/ kap); % d^-1, specific growth rate                                      
    end
    pC  = E * (vT/ L - r); % J/d, mobilized energy flux
    % generate derivatives
    dE_H  = ((1 - kap) * pC - kT_J * E_H);     % J/d, change in cumulated energy invested in maturation
    dE_R  = 0; % J/d, change in reproduction buffer
  else
    % adults do not grow, and do not shrink
    r = 0; % d^-1, specific growth rate                                      
    pC  = E * vT/ L; % J/d, mobilized energy flux
    dE_H  = 0;     % J/d, change in cumulated energy invested in maturation
    dE_R  = pC - pT_M * L^3 - kT_J * E_Hp; % J/d, change in reproduction buffer
  end
  % generate derivatives
  dE    = pA - pC;  % J/d, change in energy in reserve
  dL    = r * L / 3;    % cm^3/d, change in structural volume

  % pack derivatives
  dELHR = [dE; dL; dE_H; dE_R]; 
end
