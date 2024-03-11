% intialization has components
%1: reading the parameters from the matlab file. Depneding of what I prefer
%I can read it from the pars_init or.mat
load('results_Tachidius_discipes.mat','par') %this is to load the parameters from .mat
%into a structure called pars
% unpack par, data, auxData

cPar = parscomp_st(par);  vars_pull(par); vars_pull(cPar);
%2: a place where you specify food and temp.
T=297.15; %temperature = 24C
f=0.9;
%3: define the time you want to model
t_max=20;
t= linspace(1,t_max,1000);
%second part: think about what i want to predict
%i.e: length. then I build a code for predicting lengtht for the time I specified (I can use code from my predict file).

pars_T = T_A;
TC = tempcorr(T, T_ref, pars_T);

pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, ~, info] = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am; % J, initial energy in egg


p = [p_Am; v; p_M; k_J; kap; kap_G; E_G; E_Hb; E_Hp];
% ELHR0 = [E_0, 1E-5, 0, 0];
pars_tp = [g, k, l_T, v_Hb, v_Hp];
[tau_p, tau_b, l_p, l_b, info] = get_tp (pars_tp, f);
L_b = L_m * l_b;  
E_b = f * E_m * L_b^3; % J, energy in reserve at birth
ELHR0 = [E_b; L_b; E_Hb; 0]; % state variables at birth
[~, ELHR] = ode45(@dget_ELHR_sbp, t, ELHR0,[], p, TC, f);
L_phys = ELHR(:,2)/ del_M; 

%third: print out results. and save results into files
data_tL = readmatrix('tL24.txt');
[t_data, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t_data,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t_data,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL24 = [t_data, Lmean]; % d, cm- time, mean length
units.tL24 = {'d', 'cm'}; label.tL24 = {'time', 'area^{(1/2)}'};
temp.tL24 = C2K(24); units.temp.tL24 = 'K'; label.temp.tL24 = 'temperature';
bibkey.tL24 = {'Vigl2023'};
stdev.tL24 =Lsd ; units.stdev.tL24 = 'cm'; label.stdev.tL24 = 'standard deviation';

plot(t, L_phys);
hold on
scatter(t_data,Lmean);
xlabel('time, d'); ylabel('length, cm'); title('24 deg C')
% fldnm = 'tL24';
% data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
% figure(10); hold on
% plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
% for i = 1:size(stdev2plot,1)
% plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
% end
% plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
% xlabel('time, d'); ylabel('length, cm'); title('24 deg C')
%               set(gca,'Fontsize',10); 
%               set(gcf,'PaperPositionMode','manual');
%               set(gcf,'PaperUnits','points'); 
%               set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
%  print('results_Tachidius_discipes_06.png', '-dpng')

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
