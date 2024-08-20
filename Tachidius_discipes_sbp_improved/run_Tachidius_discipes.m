
close all;
global pets

pets = {'Tachidius_discipes'};
% check_my_pet(pets);

estim_options('default');
estim_options('max_step_number', 2e3);
estim_options('max_fun_evals',5e3);

estim_options('pars_init_method',2);
estim_options('results_output', 5);
estim_options('method', 'nm');

estim_pars;
% mat2pars_init;

%% put here personalized figures

% code below to get the prediction given parameters in pars_init - so only
% run when parameters are in pars init. you could choose to load the
% parameters from the results file
[data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes;
[par, metaPar, txtPar] = pars_init_Tachidius_discipes(metaData);
dataNew = data;
dataNew.tL12 = linspace(data.tL12(1,1), data.tL12(end,1),100)';
dataNew.tL15 = linspace(data.tL15(1,1), data.tL15(end,1),100)';
dataNew.tL18 = linspace(data.tL18(1,1), data.tL18(end,1),100)';
dataNew.tL21 = linspace(data.tL21(1,1), data.tL21(end,1),100)';
dataNew.tL24 = linspace(data.tL24(1,1), data.tL24(end,1),100)';
[prdData, info] = predict_Tachidius_discipes(par, dataNew, auxData);

% here is the code to plot the lengths with the sd, I took same color
% coding for the colors as in the automatized figure
fldnm = 'tL12';
color = [0 0 0]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(6); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('12 deg C')
              set(gca,'Fontsize',10); 
              set(gcf,'PaperPositionMode','manual');
              set(gcf,'PaperUnits','points'); 
              set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
 print('results_Tachidius_discipes_06.png', '-dpng')

fldnm = 'tL15';
color = [0, 0, 1]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(7); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('15 deg C')
              set(gca,'Fontsize',10); 
              set(gcf,'PaperPositionMode','manual');
              set(gcf,'PaperUnits','points'); 
              set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
 print('results_Tachidius_discipes_07.png', '-dpng')


fldnm = 'tL18';
color = [1, 0, 1]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(8); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('18 deg C')
              set(gca,'Fontsize',10); 
              set(gcf,'PaperPositionMode','manual');
              set(gcf,'PaperUnits','points'); 
              set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
 print('results_Tachidius_discipes_08.png', '-dpng')

fldnm = 'tL21';
color = [1, 0, 0]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(9); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('21 deg C')
              set(gca,'Fontsize',10); 
              set(gcf,'PaperPositionMode','manual');
              set(gcf,'PaperUnits','points'); 
              set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
 print('results_Tachidius_discipes_09.png', '-dpng')


fldnm = 'tL24';
color = [1, .5, .5]; 
data2plot =  data.(fldnm); stdev2plot = auxData.stdev.(fldnm); prdData2plot = prdData.(fldnm);
figure(10); hold on
plot(data2plot(:,1), data2plot(:,2), 'color',color, 'markersize', 8,'marker', 'o')
for i = 1:size(stdev2plot,1)
plot([data2plot(i,1), data2plot(i,1)],[data2plot(i,2) + stdev2plot, data2plot(i,2) - stdev2plot], 'LineStyle','-', 'color', color, 'LineWidth',1)
end
plot(dataNew.(fldnm), prdData2plot, 'color',color, 'LineStyle', '-','LineWidth', 2)
xlabel('time, d'); ylabel('length, cm'); title('24 deg C')
              set(gca,'Fontsize',10); 
              set(gcf,'PaperPositionMode','manual');
              set(gcf,'PaperUnits','points'); 
              set(gcf,'PaperPosition',[0 0 350 250]); % left bottom width height
 print('results_Tachidius_discipes_10.png', '-dpng')


prt_report_AmPtox