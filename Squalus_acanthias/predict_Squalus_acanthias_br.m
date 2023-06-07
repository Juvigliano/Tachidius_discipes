function [prdData, info] = predict_Squalus_acanthias_br(par, data, auxData)
% file generated by prt_predict

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
v2struct(par); v2struct(cPar); v2struct(data); v2struct(auxData);

% compute temperature correction factors
pars_T = T_A;
TC_ab = tempcorr(temp.ab, T_ref, pars_T);
TC_tp = tempcorr(temp.tp, T_ref, pars_T);
TC_am = tempcorr(temp.am, T_ref, pars_T);
TC_Ri = tempcorr(temp.Ri, T_ref, pars_T);
TC_LR = tempcorr(temp.LR, T_ref, pars_T);

% life cycle
pars_tp = [g k l_T v_Hb v_Hp]; % compose par vector
[tau_p, tau_b, l_p, l_b, info] = get_tp(pars_tp, f); % -, scaled times and lengths

if info == 0
  prdData = []; return;
end

% initial
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
Ww_0 = U_E0 * p_Am * w_E/ mu_E/ d_V; % g, initial wet weight

% birth
L_b = L_m * l_b;                  % cm, structural length at birth
a_b = t_0 + tau_b/ k_M/ TC_ab;    % d, age at birth
Lw_b = L_b/ del_M;                % cm, physical length at birth

% puberty
L_p = L_m * l_p; % cm, structural length at puberty
t_p = (tau_p - tau_b)/ k_M/ TC_tp; % d, time since birth at puberty
Lw_p = L_p/ del_M; % cm, physical length at puberty
Ww_p = L_p^3 * (1 + f * ome); % g, wet weight at puberty

% ultimate
l_i = f - l_T; % -, scaled ultimate length
L_i = L_m * l_i; % cm, ultimate structural length
Lw_i = f * L_m/ del_M; % cm, ultimate physical length
Ww_i = L_i^3 * (1 + f * ome); % g, ultimate wet weight 
pars_tm = [g; l_T; h_a/ k_M^2; s_G]; % compose parameter vector
tau_m = get_tm_s(pars_tm, f, l_b); % -, scaled mean life span
a_m = tau_m/ k_M/ TC_am; % d, mean life span

% reproduction
pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector
R_i = TC_Ri * reprod_rate(L_i, f, pars_R); % #/d, ultimate reproduction rate

% pack to output
prdData.ab = a_b;
prdData.tp = t_p;
prdData.am = a_m;
prdData.Lb = Lw_b;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;
prdData.Ww0 = Ww_0;
prdData.Wwp = Ww_p;
prdData.Wwi = Ww_i;
prdData.Ri = R_i;

% length - reproduction rate
pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hp]; % compose parameter vector
LR = TC_LR * reprod_rate(LR(:,1) * del_M, f_LR, pars_R); % #/d, reproduction rate

% pack to output
prdData.LR = LR;

