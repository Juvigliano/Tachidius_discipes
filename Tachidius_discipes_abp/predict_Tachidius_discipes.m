function [prdData, info] = predict_Tachidius_discipes(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par);  vars_pull(par);
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  % compute temperature correction factors
  TC     = tempcorr(temp.ah, T_ref, T_A);  kT_M = k_M * TC; % 22 C
  TC_T   = tempcorr(C2K(Tah(:,1)), T_ref, T_A);             % 14.86 C, 19.34 C, 24.90 C
  TC_Tab = tempcorr(C2K(Tah(:,1)), T_ref, T_A);
  TC_Ttj = tempcorr(C2K(Ttj(:,1)), T_ref, T_A);
  TC_Ttp = tempcorr(C2K(Ttp(:,1)), T_ref, T_A);
  
  % zero-variate data at f = f and T = 20 C (except for respiration at 20 C)
  
  % life cycle
  % life cycle
  pars_tjj = [g k l_T v_Hb v_Hj v_Hp];
  pars_tjp = [g k l_T v_Hb v_Hp v_Hp];
  
  % gets scaled age and length at birth and metam
  [t_j, ~, t_b, l_j, ~, l_b, ~, ~, ~, info_j] = get_tj(pars_tjj, f);
   % gets scaled age and length at puberty
  [t_p, ~, ~, l_p, ~, ~, ~, ~, ~, info_p] = get_tj(pars_tjp, f);
  info = info_j == 1 && info_p == 1;
  
  % initial (used in R_i)
  pars_UE0 = [V_Hb, g, k_J, k_M, v]; % compose parameter vector
  E_0 = p_Am * initial_scaled_reserve(f, pars_UE0); % J, initial reserve
  
  % birth (start of acceleration)
  L_b  = L_m * l_b;                  % cm, structural length at birth
  Lw_b = L_b/ del_M;                 % cm, length at birth
  a_b  = t_b/ k_M;  aT_b  = a_b/ TC; % d, age at birth
  
  % metam (morphological only)
  L_j  = L_m * l_j;                  % cm, structural length at metam
  Lw_j = L_j/ del_M;                 % cm, length at metam
  aT_j = t_j/ kT_M;                  % d, age at metam
  
  % puberty (end of acceleration and growth and kappa-rule)
  L_p  = L_m * l_p;                  % cm, structural length at puberty
  Lw_p = L_p/ del_M;                 % cm, length at puberty
  aT_p = t_p/ kT_M;                  % d, age at puberty
  Ww_p = 1e6 * L_p^3 * (1 + f * w);  % mug, wet weight at puberty
  
  % ultimate
  Lw_i = Lw_p;                       % cm, ultimate length
  Wd_i = Ww_p * d_V;                 % mug, ultimate dry weight
  
   
  % pack to output
  % prdData.am  = aT_m;            % d, life span
  % prdData.ah  = (t_0 + a_b)/ TC; % d, age at hatch
  prdData.tp  = aT_p - aT_j;     % d, time since metam at puberty
  prdData.Lh  = Lw_b;            % cm, length at hatch (assumed to equal length at birth)
  prdData.Lp  = Lw_p;            % cm, length at puberty
  prdData.Li  = Lw_i;            % cm, ultimate length (same as at puberty)
  prdData.Wdi = Wd_i;            % mug, ultimate dry weight (same as at puberty)
  
  % uni-variate data
  
  % time-length at f = f_tL and T = 22 C
  [~, ~, ~, l_p, ~, l_b] = get_tj(pars_tjp, f);      % overwrite l_p and l_b for f = f_tL
  L_b  = L_m * l_b;                                     % cm, structural length at birth at f = f_tL
  L_p  = L_m * l_p;                                     % cm, structural length at puberty at f = f_tL
  rT_j = TC * v * (f/ L_b - 1/ L_m)/ (f + g);     % 1/d, specific growth rate during accelleration at f = f_tL
  ELw  = min(L_p, L_b * exp(f(:,1) * rT_j/ 3))/ del_M; % cm, length at f = f_tL


  % % temperature-development time at f = f_Tt
  % [t_j, ~, t_b] = get_tj(pars_tjj, f_Tt); % overwrite t_j and t_b for f = f_Tt
  % t_p  = get_tj(pars_tjp, f_Tt);          % overwrite t_p for f = f_Tt
  % Ea_b = (t_0 + t_b/ k_M) ./ TC_T;        % d, age at birth at all temperatures at f = f_Tt
  % Et_j = (t_j - t_b)/ k_M ./ TC_T;        % d, age at metam at all temperatures at f = f_Tt
  % Et_p = (t_p - t_j)/ k_M ./ TC_T;        % d, age at puberty at all temperatures at f = f_Tt
  % 
  %  % get scaled time and length at birth at f = 1 (parents fed ad libitum)
  % [t_b, l_b] = get_tb(pars_tjj([1 2 4]), 1); % overwrite t_b and l_b
  % 
  % birth
  aT_bX = t_b/ kT_M; % d, age at birth at f = 1 (parents fed ad libitum)
  
  % pack to output
  prdData.tL12  = ELw;           % d, time        - cm, length
  % prdData.Xtj = aT_jX - aT_bX; % mug C/mL, food - d, time since hatch at metam
  % prdData.Xtp = aT_pX - aT_jX; % mug C/mL, food - d, time since metam at puberty
  % prdData.XR  = RT_iX;         % mug C/mL, food - #/d, reprod rate
  % prdData.Tah = Ea_b;          % C, temperature - d, inter-brood time
  % prdData.Ttj = Et_j;          % C, temperature - d, time since hatch at metam
  % prdData.Ttp = Et_p;          % C, temperature - d, time since metam at puberty
  % 