
legend_crust = { ... %
{'*', 8, 1, [0 0 0], [1 1 1]}, 'Tachidius_discipes', ...
{'o', 8, 1, [0 0 0], [1 1 1]}, 'Calanoidea', ...
{'o', 8, 1, [0 0 0], [0 1 1]}, 'Crustacea', ...
{'o', 8, 1, [1 1 1], [0.8 0.8 0.8]}, 'Animalia'
};



shstat_options('default');
LiEHb = read_allStat({'L_i', 'E_Hb'});
[Li_EHb, leg] = shstat(LiEHb, legend_crust, 'Crustacea');
figure(Li_EHb)
xlabel('_{10}log ultimate struc length, L_i^\infty, cm')
ylabel('_{10}log E_H^b, J')
print -r300 -dpng Li_EHb_crustacea.png
figure(leg)
print -r300 -dpng leg_crustacea.png

