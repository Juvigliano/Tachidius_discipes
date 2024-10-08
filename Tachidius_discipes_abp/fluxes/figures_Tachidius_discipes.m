% unpack par, data, auxData
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);

cPar = parscomp_st(par);  vars_pull(par);
vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);

% Set temperature and food availability

T=C2K(24)
f =0.9006;

% predicting fluxes and other statistics
[stat, txtStat] = statistics_st('abp', par, T, f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Energy budget at birth %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fluxes: birth
fluxes = [stat.p_Sb/stat.p_Cb, stat.p_Jb/stat.p_Cb, stat.p_Gb/stat.p_Cb, stat.p_Rb/stat.p_Cb];
flux_labels = {'Somatic maintenance', 'Maturity maintenance', 'Growth', 'Maturation'}; % Corresponding labels
label_with_proportions = strcat(flux_labels, {' '}, string(round(fluxes, 2)));

% pie plot
figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
h = pie(fluxes, [1 1 0 0], label_with_proportions);
title(strcat({'Energy budget at birth at '}, string(K2C(T)), {'C and f = '}, string(f)));

% Define custom colors (RGB)
colors = [
    0.2 0.6 0.8;  % color for pS
    0.8 0.2 0.2;  % color for pJ
    0.6 0.8 0.2;  % color for pG
    0.8 0.6 0.2   % color for pR
];

% Apply custom colors to pie segments
for k = 1:2:length(h)
    h(k).FaceColor = colors((k+1)/2, :);
end


% bar plot
figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
b = bar(fluxes);
set(gca, 'XTickLabel', flux_labels);
ylabel('Proportion');
title(strcat({'Energy budget at birth at '}, string(K2C(T)), {'C and f = '}, string(f)));

% Rotate the x-axis labelst
xtickangle(45);

% Define custom colors
colors = [
    0.2 0.6 0.8;  % color for pS
    0.8 0.2 0.2;  % color for pJ
    0.6 0.8 0.2;  % color for pG
    0.8 0.6 0.2   % color for pR
];

% Apply custom colors to bars
b.FaceColor = 'flat';
b.CData = colors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Energy budget after puberty %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fluxes: ultimate
fluxes_p = [stat.p_Si/stat.p_Ci, stat.p_Ji/stat.p_Ci, stat.p_Gi/stat.p_Ci, stat.p_Ri/stat.p_Ci];
flux_labels = {'Somatic maintenance', 'Maturity maintenance', 'Growth', 'Reproduction'}; % Corresponding labels
label_with_proportions = strcat(flux_labels, {' '}, string(round(fluxes_p, 2)));

% pie plot
figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
h = pie(fluxes_p, [1 1 0 0], label_with_proportions);
title(strcat({'Energy budget after puberty at '}, string(K2C(T)), {'C and f = '}, string(f)));

% Define custom colors (RGB)
colors = [
    0.2 0.6 0.8;  % color for pS
    0.8 0.2 0.2;  % color for pJ
    0.6 0.8 0.2;  % color for pG
    0.8 0.6 0.2   % color for pR
];

% Apply custom colors to pie segments
for k = 1:2:length(h)
    h(k).FaceColor = colors((k+1)/2, :);
end

% bar plot
fig = figure('Position', [100, 100, 800, 600]); % [left, bottom, width, height]
b = bar(fluxes_p);
set(gca, 'XTickLabel', flux_labels);
ylabel('Proportion');
title(strcat({'Energy budget after puberty at '}, string(K2C(T)), {'C and f = '}, string(f)));

% Rotate the x-axis labels
xtickangle(45);

% Define custom colors
colors = [
    0.2 0.6 0.8;  % color for pS
    0.8 0.2 0.2;  % color for pJ
    0.6 0.8 0.2;  % color for pG
    0.8 0.6 0.2   % color for pR
];

% Apply custom colors to bars
b.FaceColor = 'flat';
b.CData = colors;


% to save into a table use
writetable( struct2table( stat ), 'stats_24_final_final.xlsx' );