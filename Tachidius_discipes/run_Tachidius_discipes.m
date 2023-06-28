
close all;
global pets

pets = {'Tachidius_discipes'};
check_my_pet(pets);

estim_options('default');
estim_options('max_step_number', 5e2);
estim_options('max_fun_evals',5e3);

estim_options('pars_init_method', 2);
estim_options('results_output', 3);
estim_options('method', 'no');

estim_pars;


%% put here personalized figures

% code below to get the prediction given parameters in pars_init - so only
% run when parameters are in pars init. you could choose to load the
% parameters from the results file
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);
[prdData, info] = predict_Tachidius_discipes(par, data, auxData);

% here is the code to plot the lengths with the sd, I took same color
% coding for the colors as in the automatized figure
fldnm = 'tL12';
color = [0 0 0]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(50); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(data2plot(:,1), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('12 deg C')

fldnm = 'tL15';
color = [0, 0, 1]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(51); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(data2plot(:,1), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('15 deg C')

fldnm = 'tL18';
color = [1, 0, 1]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(52); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(data2plot(:,1), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('18 deg C')

fldnm = 'tL21';
color = [1, 0, 0]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(53); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(data2plot(:,1), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('21 deg C')

fldnm = 'tL24';
color = [1, .5, .5]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(54); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(data2plot(:,1), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('24 deg C')


