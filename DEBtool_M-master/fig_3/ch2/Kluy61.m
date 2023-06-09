%% fig:Kluy61
%% bib:Kluy61,Smit91,Grun87
%% out:kluy61,grun87

%% weight and feeding at age for Parus atricapillus

%% weight^1/3 (g^1/3) at age (d): \cite{Kluy61,Smit91}
aL = [1 1.118688942;
      3 1.488805553;
      5 1.754410643;
      7 1.903778262;
      9 2.064560231;
      11 2.168702885;
      13 2.223980091];

%% feeding rate (g/d) at age (d): female and male: \cite{Grund87,Smit91}
aF_f = [1 0.1055627246 0.4525488819;
      2 0.1438670505 0.7270538913;
      3 0.3164787408 0.9788878199;
      4 0.7575524603 1.419403051;
      5 1.131647025 1.69950602;
      6 1.682528864 1.898082986;
      7 2.114850051 2.687539084;
      8 2.591967989 3.154438457;
      9 2.162770797 3.43909222;
      10 2.950341996 3.396847898;
      11 2.677145895 3.517598847;
      12 3.166786161 4.595325802;
      13 3.549842691 4.636505457;
      14 3.041702505 4.112548748;
      15 3.750400521 4.391185685;
      16 3.725368943 4.46630215;
      17 3.326687821 3.699499633;
      18 3.548424546 3.370450484;
      19 2.92914983 3.202092842;
      20 2.803272099 3.533498064;
	21 2.882729925 1.523961879];


aF_m = [21 1.52397;
	20 3.53350;
	19 3.20210;
	18 3.37045;
	17 3.69951;
	16 4.46630;
	15 4.39119;
	14 4.11255;
	13 4.63651;
	12 4.59534;
	11 3.51759;
	10 3.39686;
	 9 3.43909;
	 8 3.15445;
	 7 2.68755;
	 6 1.89808;
	 5 1.69951;
	 4 1.41941;
	 3 0.97889;
	 2 0.72706;
	 1 0.45255];

nrregr_options('report',0)
pbert = nrregr('bert',[.8 2.4 .16]', aL);
pfeeding_f =  nrregr('fbert',[0 0; 2 1; .3 1], aF_f);
pfeeding_m =  nrregr('fbert',[0 0; 2 1; .4 1], aF_m);
nrregr_options('report',1)

%% gset term postscript color solid 'Times-Roman' 35

subplot(1,2,1); 
a = linspace(0,14,100)';
plot(aL(:,1), aL(:,2),'.g', a, bert(pbert, a),'-r')
xlabel('age, d');
ylabel('weight^1/3, g^1/3');
		    
subplot(1,2,2); 
aF = linspace(0,21,100)';
F_f = fbert (pfeeding_f(:,1), aF);
F_m = fbert (pfeeding_m(:,1), aF);
plot (aF, F_f, '-r', aF, F_m, '-b',  ...
      aF_m(:,1), aF_m(:,2), '.b' , aF_f(:,1), aF_f(:,2), '.r');
legend('female', 'male');
xlabel('age, d');
ylabel('feeding rate, g/d');
