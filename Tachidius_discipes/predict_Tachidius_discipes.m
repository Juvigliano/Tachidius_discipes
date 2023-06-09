function [prdData, info] = predict_Tachidius_discipes(par, data, auxData)
% file generated by prt_predict

% unpack par, data, auxData
cPar = parscomp_st(par); vars_pull(par);
v2struct(par); v2struct(cPar); v2struct(data); v2struct(auxData);

filterChecks= E_Hh>= E_Hb;

if filterChecks
    info=0;
    prdData={};
    return;
end

% compute temperature correction factors
pars_T = T_A;
TC_ah = tempcorr(temp.ah, T_ref, pars_T);
TC_tL = tempcorr(temp.tL, T_ref, pars_T);
TC_tL_18 = tempcorr(temp.tL_18, T_ref, pars_T);
TC_tL_12 = tempcorr(temp.tL_12, T_ref, pars_T);
TC_tL_21 = tempcorr(temp.tL_21, T_ref, pars_T);
TC_tL_24 = tempcorr(temp.tL_24, T_ref, pars_T);

TC_tN_12 = tempcorr(temp.tN_12, T_ref, pars_T);
TC_tN_15 = tempcorr(temp.tN_15, T_ref, pars_T);
TC_tN_18 = tempcorr(temp.tN_18, T_ref, pars_T);
TC_tN_21 = tempcorr(temp.tL_21, T_ref, pars_T);
TC_tN_24 = tempcorr(temp.tL_24, T_ref, pars_T);

% life cycle
pars_tj = [g; k; l_T; v_Hb; v_Hp-1e-6; v_Hp];
[tau_j, tau_p, tau_b, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_tj(pars_tj, f);

if info == 0
  prdData = []; return;
end


% hatch
pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
U_E0 = initial_scaled_reserve(f, pars_UE0); % d.cm^2, initial scaled reserve
[U_H aUL] = ode45(@dget_aul, [0; U_Hh], [0 U_E0 1e-10], [], kap, v, k_J, g, L_m);
a_h = aUL(end,1)/ TC_ah; % d, age at hatch
Lw_h = aUL(end,3)/ del_M; % cm, physical length at hatch

% birth
L_b = L_m * l_b;                  % cm, structural length at birth
Lw_b = L_b/ del_M;                % cm, physical length at birth

% puberty
L_p = L_m * l_p; % cm, structural length at puberty
Lw_p = L_p/ del_M; % cm, physical length at puberty

% ultimate
l_i = f - l_T; % -, scaled ultimate length
L_i = L_m * l_i; % cm, ultimate structural length
Lw_i = L_p/ del_M; % cm, ultimate physical length
Wd_i = L_p^3 * (1 + f * ome) * d_V; % g, ultimate dry weight

% pack to output
prdData.ah = a_h;
prdData.Lh = Lw_h;
prdData.Lb = Lw_b;
prdData.Lp = Lw_p;
prdData.Li = Lw_i;
prdData.Wdi = Wd_i;

% time - length
L_b = L_m * get_lb([g, k, v_Hb], f_tL); L_i = L_m * f_tL; % cm, structural lengths
rT_B = TC_tL * k_M/ 3/ (1 + f_tL/ g); % 1/d, von Bert growth rate
EL = (L_i - (L_i - L_b) * exp( - rT_B * tL(:,1)))/ del_M; % cm, physical length

% time - length
L_b = L_m * get_lb([g, k, v_Hb], f_tL); L_i = L_m * f_tL; % cm, structural lengths
rT_B = TC_tL_18 * k_M/ 3/ (1 + f_tL/ g); % 1/d, von Bert growth rate
EL18 = (L_i - (L_i - L_b) * exp( - rT_B * tL_18(:,1)))/ del_M; % cm, physical length

% time - length
L_b = L_m * get_lb([g, k, v_Hb], f_tL); L_i = L_m * f_tL; % cm, structural lengths
rT_B = TC_tL_12 * k_M/ 3/ (1 + f_tL/ g); % 1/d, von Bert growth rate
EL12 = (L_i - (L_i - L_b) * exp( - rT_B * tL_12(:,1)))/ del_M; % cm, physical length

% time - length
L_b = L_m * get_lb([g, k, v_Hb], f_tL); L_i = L_m * f_tL; % cm, structural lengths
rT_B = TC_tL_21 * k_M/ 3/ (1 + f_tL/ g); % 1/d, von Bert growth rate
EL21 = (L_i - (L_i - L_b) * exp( - rT_B * tL_21(:,1)))/ del_M; % cm, physical length

% time - length
L_b = L_m * get_lb([g, k, v_Hb], f_tL); L_i = L_m * f_tL; % cm, structural lengths
rT_B = TC_tL_24 * k_M/ 3/ (1 + f_tL/ g); % 1/d, von Bert growth rate
EL24 = (L_i - (L_i - L_b) * exp( - rT_B * tL_24(:,1)))/ del_M; % cm, physical length

 % tN data 12 C
  pars_R = [kap; kap_R; g; k_J*TC_tN_12; k_M*TC_tN_12; L_T; v*TC_tN_12; U_Hb/TC_tN_12; U_Hp/TC_tN_12]; % pars for cum_reprod
EN = cum_reprod([0; tN_12(:,1)], f, pars_R);
EN_12=ones(length(tN_12(:,1)),1)*EN(end);

 % tN data 15 C
  pars_R = [kap; kap_R; g; k_J*TC_tN_15; k_M*TC_tN_15; L_T; v*TC_tN_15; U_Hb/TC_tN_15; U_Hp/TC_tN_15]; % pars for cum_reprod
EN = cum_reprod([0; tN_15(:,1)], f, pars_R);
EN_15=ones(length(tN_15(:,1)),1)*EN(end);


 % tN data 18 C
  pars_R = [kap; kap_R; g; k_J*TC_tN_18; k_M*TC_tN_18; L_T; v*TC_tN_18; U_Hb/TC_tN_18; U_Hp/TC_tN_18]; % pars for cum_reprod
EN = cum_reprod([0; tN_18(:,1)], f, pars_R);
EN_18=ones(length(tN_18(:,1)),1)*EN(end);

% tN data 21 C
  pars_R = [kap; kap_R; g; k_J*TC_tN_21; k_M*TC_tN_21; L_T; v*TC_tN_21; U_Hb/TC_tN_21; U_Hp/TC_tN_21]; % pars for cum_reprod
EN = cum_reprod([0; tN_21(:,1)], f, pars_R);
EN_21=ones(length(tN_21(:,1)),1)*EN(end);


% tN data 24 C
  pars_R = [kap; kap_R; g; k_J*TC_tN_24; k_M*TC_tN_24; L_T; v*TC_tN_24; U_Hb/TC_tN_24; U_Hp/TC_tN_24]; % pars for cum_reprod
EN = cum_reprod([0; tN_24(:,1)], f, pars_R);
EN_24=ones(length(tN_24(:,1)),1)*EN(end);


% pack to output
prdData.tL = EL;
prdData.tL_18 = EL18;
prdData.tL_12 = EL12;
prdData.tL_21 = EL21;
prdData.tL_24 = EL24;

% pack to output
prdData.tN_12= EN_12;
prdData.tN_15 = EN_15;
prdData.tN_18 = EN_18;
prdData.tN_21 = EN_21;
prdData.tN_24 = EN_24;