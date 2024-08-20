% intialization has components
%1: reading the parameters from the matlab file. Depneding of what I prefer
%I can read it from the pars_init or.mat
load('results_Tachidius_discipes.mat','par') %this is to load the parameters from .mat
%into a structure called pars
% unpack par, data, auxData

cPar = parscomp_st(par);  vars_pull(par); vars_pull(cPar);
%2: a place where you specify food and temp.
T=294.15; %temperature = 21C
%3: define the time you want to model
t_max=30;
t= linspace(1,t_max,1000);
%second part: think about what i want to predict
%i.e: length. then I build a code for predicting lengtht for the time I specified (I can use code from my predict file).

pars_T = T_A;
TC = tempcorr(T, T_ref, pars_T);

pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
[U_E0, ~, info] = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
E_0 = U_E0 * p_Am; % J, initial energy in egg

% Define a range of f values
f_values = [0.85, 0.9, 1];

% Initialize a cell array to store model predictions for different f values
L_phys_all = cell(length(f_values), 1);

% Loop over each f value
for i = 1:length(f_values)
    f = f_values(i);
    
    % Re-run the model to obtain predictions of length over time
    [~, ELHR] = ode45(@dget_ELHR_sbp, t, ELHR0,[], p, TC, f);
    L_phys_all{i} = ELHR(:,2) / del_M; % Store the model predictions for the current f value
end
%third: print out results. and save results into files
data_tL = readmatrix('tL21.txt');
[t_data, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t_data,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t_data,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL21 = [t_data, Lmean]; % d, cm- time, mean length
units.tL21 = {'d', 'cm'}; label.tL21 = {'time', 'area^{(1/2)}'};
temp.tL21 = C2K(21); units.temp.tL21 = 'K'; label.temp.tL21 = 'temperature';
bibkey.tL21 = {'Vigl2023'};
stdev.tL21 =Lsd ; units.stdev.tL21 = 'cm'; label.stdev.tL21 = 'standard deviation';
% Plot model predictions for different f values
figure;
hold on;
for i = 1:length(f_values)
    plot(t, L_phys_all{i}, 'DisplayName', ['f = ', num2str(f_values(i))]);
end
xlabel('Time (days)');
ylabel('Length (cm)');
title('Model Predictions for Different f Values at 21  C');
legend('show');
hold on
scatter(t_data,Lmean);
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