
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


  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  [U_E0, ~, info] = initial_scaled_reserve(1.0, pars_UE0); % d.cm^2, initial scaled reserve
  if info == 0;  prdData = []; return; end
  E_0 = U_E0 * p_Am; % J, initial energy in egg

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

%%time- number of offspring
% data modified to work with means per day
data.tN12 = [...
   37.00	mean([25,29,22,26,33,33,27,25,30,29,20])];
units.tN12 = {'d', '#'}; label.tN12= {'time', 'clutch size'};
temp.tN12 = C2K(12); units.temp.tN12 = 'K'; label.temp.tN12 = 'temperature';
bibkey.tN12 = {'Vigl2023'}; treat.tN12 = {0};
%
data.tN15 = [...
32.00	mean([41,25,10])
33.00	31
37.00	mean([33,31,21,25])];
units.tN15 = {'d', '#'}; label.tN15= {'time', 'clutch size'};
temp.tN15 = C2K(15); units.temp.tN15 = 'K'; label.temp.tN15 = 'temperature';
bibkey.tN15 = {'Vigl2023'};
treat.tN15 = {0};
%
data.tN18 = [...
23.00	18
24.00	mean([20,18,29, 21])
30.00	mean([24,34,34,32,46])];
units.tN18 = {'d', '#'}; label.tN18= {'time', 'clutch size'};
temp.tN18 = C2K(18); units.temp.tN18 = 'C'; label.temp.tN18 = 'temperature';
bibkey.tN18 = {'Vigl2023'};
treat.tN18 = {0};
%
data.tN21 = [...
17.00	29
18.00	18
20.00	5
23.00	18
22.00	mean([17,12])
23.00	mean([24,34])
24.00	10];
units.tN21 = {'d', '#'}; label.tN21= {'time', 'clutch size'};
temp.tN21 = C2K(21); units.temp.tN21 = 'K'; label.temp.tN21 = 'temperature';
bibkey.tN21 = {'Vigl2023'};
treat.tN21 = {0};
%
data.tN24 = [...
17.00	12
18.00	21
19.00	mean([10,23,10])
18.00	mean([19,11,12,	11,	10,	16,	15])];
units.tN24 = {'d', '#'}; label.tN24= {'time', 'clutch size'};
temp.tN24 = C2K(24); units.temp.tN24 = 'C'; label.temp.tN24 = 'temperature';
bibkey.tN24 = {'Vigl2023'};
treat.tN24 = {0}; 

% time-clutch size at 12 C
  EN12 = zeros(size(tN12,1),1);
  for i = 1:size(tN12,1)
      [~, ELHR] = ode45(@dget_ELHR_abp, [0 tN12(i,1)], ELHR0,[], p, TC12, L_b, f_tL12);
      EN12(i) = kap_R * ELHR(end,4)/E_0; 
  end

  EN15 = zeros(size(tN15,1),1);
  for i = 1:size(tN15,1)
      [~, ELHR] = ode45(@dget_ELHR_abp, [0 tN15(i,1)], ELHR0,[], p, TC15, L_b, f_tL15);
      EN15(i) = kap_R * ELHR(end,4)/E_0; 
  end

  EN18 = zeros(size(tN18,1),1);
  for i = 1:size(tN18,1)
      [~, ELHR] = ode45(@dget_ELHR_abp, [0 tN18(i,1)], ELHR0,[], p, TC18, L_b, f_tL18);
      EN18(i) = kap_R * ELHR(end,4)/E_0; 
  end

  EN21 = zeros(size(tN21,1),1);
  for i = 1:size(tN21,1)
      [~, ELHR] = ode45(@dget_ELHR_abp, [0 tN21(i,1)], ELHR0,[], p, TC21, L_b, f_tL21);
      EN21(i) = kap_R * ELHR(end,4)/E_0; 
  end

  EN24 = zeros(size(tN24,1),1);
  for i = 1:size(tN24,1)
      [~, ELHR] = ode45(@dget_ELHR_abp, [0 tN24(i,1)], ELHR0,[], p, TC24, L_b, f_tL21);
      EN24(i) = kap_R * ELHR(end,4)/E_0; 
  end
figure;
hold on;

% Extract time and clutch size for each temperature
tN12_time = data.tN12(:,1);
tN12_size = data.tN12(:,2);
tN15_time = data.tN15(:,1);
tN15_size = data.tN15(:,2);
tN18_time = data.tN18(:,1);
tN18_size = data.tN18(:,2);
tN21_time = data.tN21(:,1);
tN21_size = data.tN21(:,2);
tN24_time = data.tN24(:,1);
tN24_size = data.tN24(:,2);

% Define colors for each temperature
colors = lines(5); % Generates a colormap with 5 distinct colors

% Scatter plot with observed values
% Define the color
color12=("#4575b4")
 color15=("#74add1")
 color18=("#abd9e9")
 color21=("#fc8d59")
 color24=("#d73027")

% Plot predictions for each temperature
plot(tN12_time, EN12,'Color', color12, 'LineWidth', 2, 'DisplayName', '12°C Prediction');scatter(tN12_time, tN12_size, 50, 'MarkerEdgeColor',color12);
plot(tN15_time, EN15, 'Color', color15,  'LineWidth', 2, 'DisplayName', '15°C Prediction');scatter(tN15_time, tN15_size, 50, 'MarkerEdgeColor',color15);
plot(tN18_time, EN18, 'Color', color18, 'LineWidth', 2, 'DisplayName', '18°C Prediction');scatter(tN18_time, tN18_size, 50, 'MarkerEdgeColor',color18);
plot(tN21_time, EN21, 'Color', color21, 'LineWidth', 2, 'DisplayName', '21°C Prediction');scatter(tN21_time, tN21_size, 50, 'MarkerEdgeColor',color21);
plot(tN24_time, EN24, 'Color', color24, 'LineWidth', 2, 'DisplayName', '24°C Prediction');scatter(tN24_time, tN24_size, 50, 'MarkerEdgeColor',color24);

% Add labels, title, and legend
xlabel('Time (days)','FontSize', 15);
ylabel('Clutch Size (#)','FontSize', 15);
% title('Time vs. Clutch Size at Different Temperatures');
legend('show'); % Automatically display the legend
grid on;
hold off;

% Create tables for each temperature
table12 = table(tN12_time, EN12, 'VariableNames', {'Time', 'Prediction'});
table15 = table(tN15_time, EN15, 'VariableNames', {'Time', 'Prediction'});
table18 = table(tN18_time, EN18, 'VariableNames', {'Time', 'Prediction'});
table21 = table(tN21_time, EN21, 'VariableNames', {'Time', 'Prediction'});
table24 = table(tN24_time, EN24, 'VariableNames', {'Time', 'Prediction'});

% Optionally, save the table to a file
writetable(table12, 'clutch_size12.csv');
writetable(table15, 'clutch_size15.csv');
writetable(table18, 'clutch_size18.csv');
writetable(table21, 'clutch_size21.csv');
writetable(table24, 'clutch_size24.csv');
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