
% unpack par, data, auxData
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);

cPar = parscomp_st(par);  vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);


%2: a place where you specify food and temp.
pars_T = [T_A, T_L, T_H, T_AL, T_AH]; % with all tolerance limits

TC12 = tempcorr(temp.tL12, T_ref, pars_T);                % temp corr for both length and clutch size
  TC15 = tempcorr(temp.tL15, T_ref, pars_T);
  TC18 = tempcorr(temp.tL18, T_ref, pars_T);
  TC21 = tempcorr(temp.tL21, T_ref, pars_T);
  TC24 = tempcorr(temp.tL24, T_ref, pars_T);
% T=285.15; %temperature = 12C
%3: define the time you want to model
t_max=30;
t= linspace(1,t_max,1000);
%second part: think about what i want to predict
%i.e: length. then I build a code for predicting lengtht for the time I specified (I can use code from my predict file).


% birth (start of acceleration)
% life cycle (modified from previous version which was pars_tp)
  pars_tjp = [g k l_T v_Hb v_Hp v_Hp];
 [~, ~, ~, l_p, ~, l_b, ~, ~, ~, info] = get_tj(pars_tjp, f);
  L_b  = L_m * l_b;                  % cm, structural length at birth


 p = [p_Am; v; p_M; k_J; kap; kap_G; E_G; E_Hb; E_Hp];
 E_b = E_m * L_b^3; % J, energy in reserve at birth
 ELHR0 = [E_b; L_b; E_Hb; 0]; % state variables at birth



% Define a range of f values
par.f_tL12 = 1.1376;    
par.f_tL15 = 1;      
par.f_tL18 = 0.95834;      
par.f_tL21 = 0.90059; 


% Initialize a cell array to store model predictions for different f values

data_tL = readmatrix('tL12.txt');
[t_data, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t_data,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t_data,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end


data_tL15 = readmatrix('tL15.txt');
[t_data15, ai, ci] = unique(data_tL15 (:,1)); maxci = max(ci); 
Lmean_15 = zeros(size(t_data15,1),1); % preallocate zeros for mean L values
Lsd_15 = zeros(size(t_data15,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean_15(i) =    mean(data_tL15 (ci==i,2));
 Lsd_15(i) =    std(data_tL15 (ci==i,2));
end

data_tL18 = readmatrix('tL18.txt');
[t_data18, ai, ci] = unique(data_tL18 (:,1)); maxci = max(ci); 
Lmean_18 = zeros(size(t_data18,1),1); % preallocate zeros for mean L values
Lsd_18 = zeros(size(t_data18,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean_18(i) =    mean(data_tL18 (ci==i,2));
 Lsd_18(i) =    std(data_tL18 (ci==i,2));
end
 
data_tL21 = readmatrix('tL21.txt');
[t_data21, ai, ci] = unique(data_tL21 (:,1)); maxci = max(ci); 
Lmean_21 = zeros(size(t_data21,1),1); % preallocate zeros for mean L values
Lsd_21 = zeros(size(t_data21,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean_21(i) =    mean(data_tL21 (ci==i,2));
 Lsd_21(i) =    std(data_tL21 (ci==i,2));
end
 
data_tL24 = readmatrix('tL24.txt');
[t_data24, ai, ci] = unique(data_tL24 (:,1)); maxci = max(ci); 
Lmean_24 = zeros(size(t_data24,1),1); % preallocate zeros for mean L values
Lsd_24 = zeros(size(t_data24,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean_24(i) =    mean(data_tL24 (ci==i,2));
 Lsd_24(i) =    std(data_tL24(ci==i,2));
end


% solving length-time
[tsort, ~, ci] = unique(tL12(:,1)); 
  [~, ELHR] = ode45(@dget_ELHR_abp, tsort, ELHR0,[], p, TC12, L_b, f_tL12);
  EL12 = ELHR(ci,2)./ del_M; 


[tsort, ~, ci] = unique(tL15(:,1)); 
  [~, ELHR] = ode45(@dget_ELHR_abp, tsort, ELHR0,[], p, TC15, L_b, f_tL15);
  EL15 = ELHR(ci,2)/ del_M; 

  [tsort, ~, ci] = unique(tL18(:,1)); 
  [~, ELHR] = ode45(@dget_ELHR_abp, tsort, ELHR0,[], p, TC18, L_b, f_tL18);
  EL18 = ELHR(ci,2)/ del_M; 

  [tsort, ~, ci] = unique(tL21(:,1)); 
  [~, ELHR] = ode45(@dget_ELHR_abp, tsort, ELHR0,[], p, TC21, L_b, f_tL21);
  EL21 = ELHR(ci,2)/ del_M; 

  [tsort, ~, ci] = unique(tL24(:,1)); 
  [~, ELHR] = ode45(@dget_ELHR_abp, tsort, ELHR0,[], p, TC24, L_b, f_tL21);
  EL24 = ELHR(ci,2)/ del_M; 

figure;
hold on;
% Define the color
color12=("#4575b4")
 color15=("#74add1")
 color18=("#abd9e9")
 color21=("#fc8d59")
 color24=("#d73027")

    plot(EL12,'Color', color12, 'LineWidth', 2,'DisplayName', '12°C Prediction'); scatter(t_data,Lmean,'MarkerEdgeColor',color12 );
    plot(EL15,'Color', color15, 'LineWidth',2, 'DisplayName', '15°C Prediction'); scatter(t_data15,Lmean_15,'MarkerEdgeColor',color15);
    plot(EL18,'Color', color18, 'LineWidth',2,'DisplayName', '18°C Prediction'); scatter(t_data18,Lmean_18,'MarkerEdgeColor',color18);
    plot(EL21,'Color', color21, 'LineWidth',2,'DisplayName', '21°C Prediction'); scatter(t_data21,Lmean_21,'MarkerEdgeColor',color21);
    plot(EL24,'Color', color24, 'LineWidth',2, 'DisplayName', '24°C Prediction'); scatter(t_data24,Lmean_24,'MarkerEdgeColor',color24);
    % 

% Add labels, title, and legend
xlabel('Age (days)','FontSize', 15);
ylabel('(Top view area)^{(1/2)}','FontSize', 15);
legend('show'); % Automatically display the legend
grid on;
hold off;
% Save the figure as a JPEG file
saveas(gcf, 'length_time.jpeg');

function dELHR = dget_ELHR_abp(t, ELHR, p, TC, L_b, f)
  % Define changes in the state variables for abp model
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

  % temp and acceleration correction
  if E_H < E_Hb
      s_M = 1;
  else
      s_M = L / L_b;
  end
  pT_Am = TC * p_Am * s_M;
  vT = TC * v * s_M;  
  pT_M = TC * p_M;
  kT_J = TC * k_J; 

  pA = (pT_Am * f * L^2) * (E_H >= E_Hb);
  if E_H < E_Hp
    if  kap * E * vT >= pT_M * L^4 % section 4.1.5 comments to Kooy2010 (if the energy is higher than mantainance costs)
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