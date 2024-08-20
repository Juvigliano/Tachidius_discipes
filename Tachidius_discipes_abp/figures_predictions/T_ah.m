
% unpack par, data, auxData
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);

cPar = parscomp_st(par);  vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);


%2: a place where you specify food and temp.
pars_T = [T_A, T_L, T_H, T_AL, T_AH]; % with all tolerance limits

TC_Tah = tempcorr(C2K(Tah(:,1)), T_ref, pars_T); % 12, 15, 18, 21, 24 C

% temperature-age at hatch
data.Tah= [...
24	2
24	3
24	4
24	4
24	3
24	4
24	2
24	2
24	2
24	2
24	2
24	2
21	3
21	2
21	1
21	3
21	2
21	2
21	2
21	1
21	1
18	2
18	3
18	3
18	3
18	1
18	4
18	4
18	3
18	4
18	4
15	5
15	5
15	3
15	4
15	5
15	5
15	4
15	4
12	7
12	7
12	5
12	5
12	5
12	5
12	5
12	5
12	5
12	5
12	5
12	7
];
units.Tah = {'K', 'd'}; label.Tah= {'temperature', 'age at hatching'};
bibkey.Tah = {'Vigl2023'};
% 

 % Prediction temperature-age at hatch
  [~, aUL] = ode45(@dget_aul, [0; U_Hh], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
  ETah = aUL(end,1)./ TC_Tah; % d, age at hatch
 
  %output
  prdData.Tah   = ETah;           % C, temperature - d, age at hatch
  
    % Define the color
  color12=("#4575b4")
figure;
hold on;
  %plot
  plot(Tah(:,1),ETah,'Color', color12,'LineWidth', 2, 'DisplayName', 'Prediction'); 
  scatter(Tah(:,1), Tah(:,2), 'MarkerEdgeColor',color12);


% Add labels, title, and legend
xlabel('Temperature (°C)','FontSize', 15);
ylabel('Age at hatching (days)','FontSize', 15);
% title('Time vs. Clutch Size at Different Temperatures');
legend('show'); % Automatically display the legend
grid on;
hold off;


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