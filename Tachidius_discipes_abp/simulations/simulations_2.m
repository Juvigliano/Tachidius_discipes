
% % Code to do simulations at the different temp scenarios and at the
% % estimated f
%load('results_Tachidius_discipes.mat'); simu_my_pet({metaData, metaPar, par}, [0 C2K(12); 30 C2K(12)] ,1.13)
% ;
% 
% load('results_Tachidius_discipes.mat'); simu_my_pet({metaData, metaPar, par}, [0 C2K(15); 30 C2K(15)] ,1)
% ;
% 
% load('results_Tachidius_discipes.mat'); simu_my_pet({metaData, metaPar, par}, [0 C2K(18); 30 C2K(18)] ,0.95)
% ;
% 
% load('results_Tachidius_discipes.mat'); simu_my_pet({metaData, metaPar, par}, [0 C2K(21); 30 C2K(21)] ,0.90)
% ;
% load('results_Tachidius_discipes.mat'); simu_my_pet({metaData, metaPar, par}, [0 C2K(24); 30 C2K(24)] ,0.9)
% ;


% Putting a heat wave in the middle with my own data
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);


% If I want to model sin functions of f and T
% t = linspace(0, 30, 30); tT_0 = 90; tf_0 =0;
% T = C2K(24) + 10* sin(1 * pi * (t + tT_0)/ 10);
% f = 0.8 + 0.2 * sin(1 * pi * (t + tf_0)/10);
% tT =[t',T']; tf =[t',f'];


% If I want to import my own tT data
tT =[...
1	10
2	10
3	10
4	10
5	10
6	10
7	10
8	10
9	10
10	10
11	11
12	12
13	13
14	14
15	16
16	18
17	20
18	24
19	24
20	24
21	25
22	24
23	24
24	23
25	22
26	21
27	18
28	17
29	16
30	11];
units.tT = {'d', 'C'}; label.tT= {'time', 'Temperature'};

t = tT(:, 1); % Time (days)
T = C2K(tT(:, 2)); % Temperature (Kelvin)
tT = [t, T]; % Combine time and temperature

t = (1:30)';
f=1;
tf =[t,f * ones(size(t))];

model = metaPar.model;
switch model
    case {'stf', 'stx', 'hep', 'hex', 'hax'}
    fprintf('Warning from get_inDyn_mod: %s model is not coded yet; work in progress\n', model);
    return
end


% temperature
if ~exist('tT','var') || isempty(tT) 
  tT = [0 metaData.T_typical];
elseif size(tT,2) == 1
  tT = [0 tT(1)];
end

% supply food 
if ~exist('tf','var') || isempty(tf)
  tf = [0 1];  
elseif size(tf,2) == 1
  tf = [0 tf(1)];
end

if size(tT,1) == 1 && size(tf,1) > 1 
    tT = [0 tT(1,2); tf(end,1) tT(1,2)];
elseif size(tT,1) > 1 && size(tf,1) == 1 
    tf = [0 tf(1,2); tT(end,1) tf(1,2)];
elseif size(tT,1) > 1 && size(tf,1) > 1  
    tSim = min(tT(end,1), tf(end,1));
    tT = tT(tT(:,1)<tSim,:); tT = [tT; tSim, tT(end,2)];
    tf = tf(tf(:,1)<tSim,:); tf = [tf; tSim, tf(end,2)];
end

[tELHR, tWNXO, tpAMGRD, aLW, aLWc]  = get_indDyn_mod(model, par, tT, tf);
n_col = size(tpAMGRD,2);


% Added from Dina
save( 'var_out.mat', 'tELHR', 'tWNXO', 'tpAMGRD', 'aLW', 'aLWc')
% writematrix([tELHR, tWNXO, tpAMGRD], 'var_out.csv')

%% plotting

close all
%
if size(tT,1) == 1 
    tT =[tELHR(:,1), ones(length(tELHR(:,1)),1)*tT(:,2)];
end
if size(tf,1) == 1 
    tf =[tELHR(:,1), ones(length(tELHR(:,1)),1)*tf(:,2)];
end

figure(1) % environment
subplot(2,1,1), hold on
plot(tT(:,1), K2C(tT(:,2)), 'k', 'Linewidth', 2)
plot([tT(1,1), tT(end,1)], [K2C(mean(tT(:,2))), K2C(mean(tT(:,2)))], 'r:', 'Linewidth', 2)
xlabel('time, d'), ylabel('temperature, C')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,1,2), hold on
plot(tf(:,1), tf(:,2), 'k', 'Linewidth', 2)
plot([tf(1,1), tf(end,1)], [mean(tf(:,2)), mean(tf(:,2))], 'r:', 'Linewidth', 2)
xlabel('time, d'), ylabel('functional response, -')
set(gca, 'FontSize', 15, 'Box', 'on')

figure(2) % t-ELHR
subplot(2,2,1), hold on
plot(tELHR(:,1), tELHR(:,2), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('reserve energy, J')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,2), hold on
plot(tELHR(:,1), tELHR(:,3), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('structural length, cm')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,3), hold on
plot(tELHR(:,1), tELHR(:,4), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('cum. energy to maturation, J')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,4), hold on
plot(tELHR(:,1), tELHR(:,5), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('cum. energy to reprod, J')
set(gca, 'FontSize', 15, 'Box', 'on')
%
figure(3) % t-LWNJX
subplot(2,2,1), hold on
plot(tWNXO(:,1), tWNXO(:,2), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('wet weight, g')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,2), hold on
plot(tWNXO(:,1), tWNXO(:,3), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('cum. number of eggs, #')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,3), hold on
plot(tWNXO(:,1), tWNXO(:,4), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('feeding rate, g/d')
set(gca, 'FontSize', 15, 'Box', 'on')
subplot(2,2,4), hold on
plot(tWNXO(:,1), tWNXO(:,5), 'k', 'Linewidth', 2)
xlabel('age, d'), ylabel('O_2 consumption, g/d')
set(gca, 'FontSize', 15, 'Box', 'on')

figure(4) % t-powers
clr = {'b-', 'k-', 'r-', 'r:', 'k:'};
hold on
for i = 2:n_col
    plot(tpAMGRD(:,1), tpAMGRD(:,i), clr{i-1}, 'Linewidth', 2)
end
legend('assim', 'maint', 'growth', 'matur/reprod', 'dissip', 'Location', 'best')
xlabel('age, d'), ylabel('powers, J/d')
set(gca, 'FontSize', 15, 'Box', 'on')
%
switch model
    case {'std', 'sbp'}
        prt_tab({{'age at birth; a_b (d)'; 'age at puberty; a_p (d)'; 'life span; a_m (d)'; ...
        'struc length at birth; L_b (cm)'; 'struc length at puberty; L_p (cm)'; 'ultimate length; L_i (cm)';  ...
        'wet weight at birth; Ww_b (g)'; 'wet weight at puberty; Ww_p (g)'; 'ultimate wet weigh; Ww_i (g)'}, aLW, aLWc}, ...
        {'description; symbol (units)', 'values at T&f (event function)', 'values at mean T&f (DEBtool functions)'});
    case 'abj'
        prt_tab({{'age at birth; a_b (d)'; 'age at metamorphosis; a_j (d)'; 'age at puberty; a_p (d)'; 'life span; a_m (d)'; ...
        'struc length at birth; L_b (cm)'; 'struc length at metamorphosis; L_j (cm)';  ...
        'struc length at puberty; L_p (cm)'; 'ultimate length; L_i (cm)';  ...
        'wet weight at birth; Ww_b (g)'; ...
        'wet weight at metamorphosis; Ww_j (g)'; 'wet weight at puberty; Ww_p (g)'; 'ultimate wet weigh; Ww_i (g)'}, aLW, aLWc}, ...
        {'description; symbol (units)', 'values at T&f (event function)', 'values at mean T&f (DEBtool functions)'});
    case 'asj'
        prt_tab({{'age at birth; a_b (d)'; 'age at start of acceleration; a_s (d)'; 'age at metamorphosis; a_j (d)'; 'age at puberty; a_p (d)'; 'life span; a_m (d)'; ...
        'struc length at birth; L_b (cm)'; 'struc length at start of acceleration; L_s (cm)'; 'struc length at metamorphosis; L_j (cm)';  ...
        'struc length at puberty; L_p (cm)'; 'ultimate length; L_i (cm)';  ...
        'wet weight at birth; Ww_b (g)'; 'wet weight at start of acceleration; Ww_s (g)';...
        'wet weight at metamorphosis; Ww_j (g)'; 'wet weight at puberty; Ww_p (g)'; 'ultimate wet weigh; Ww_i (g)'}, aLW, aLWc}, ...
        {'description; symbol (units)', 'values at T&f (event function)', 'values at mean T&f (DEBtool functions)'});
    case 'abp'
        prt_tab({{'age at birth; a_b (d)'; 'age at puberty; a_p (d)'; 'life span; a_m (d)'; ...
        'struc length at birth; L_b (cm)'; 'struc length at puberty; L_p (cm)'; 'ultimate length; L_i (cm)';  ...
        'wet weight at birth; Ww_b (g)'; 'wet weight at puberty; Ww_p (g)'; 'ultimate wet weigh; Ww_i (g)'}, aLW, aLWc}, ...
        {'description; symbol (units)', 'values at T&f (event function)', 'values at mean T&f (DEBtool functions)'});
end
