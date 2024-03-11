function [data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes
% file generated by prt_mydata

%% set metaData
metaData.phylum     = 'Arthropoda';
metaData.class      = 'Copepoda';
metaData.order      = 'Harpacticoida';
metaData.family     = 'Tachidiidae';
metaData.species    = 'Tachidius_discipes';
metaData.species_en = 'Copepod';

% metaData.ecoCode.climate = {'Dfb','Cfb'};
% metaData.ecoCode.ecozone = {'MAb','MAn'};
% metaData.ecoCode.habitat = {'jiSm','jiM','jiMm','jiMi','jiMcb'};
% metaData.ecoCode.embryo  = {'Mbf'};
% metaData.ecoCode.migrate = {'TW'};
% metaData.ecoCode.food    = {'jiD','jiB','jiPp'};
% metaData.ecoCode.gender  = {'D'};
% metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(15); % K, body temp

metaData.data_0     = {'ah','Lh','Lb','Lp','Wdp','am'};
metaData.data_1     = {'tN_T','tR_T','tL_T','LWC','LWN','T-ah','t-ap'};

metaData.COMPLETE   = 2; % using criteria of LikaKear2011

metaData.author     = {'Julieta Vigliano Relva'};
metaData.date_subm  = [2023 6 16];
metaData.email      = {'julieta.vigliano@ugent.be'};
metaData.address    = {'UGent Ghent University'};

metaData.author_mod_1    = {'Julieta Vigliano Relva'};
metaData.date_mod_1 = [2023 6 16];
metaData.email_mod_1      = {'julieta.vigliano@ugent.be'};
metaData.address_mod_1    = {'UGent Ghent University'};

% metaData.curator    = {'Starrlight Augustine'};
% metaData.email_cur  = {'starrlight@tecnico.ulisboa.pt'};
% metaData.date_acc   = [2023 6 21];

%% set zero-variate data
data.ah = 4.375; units.ah = 'd'; label.ah = 'age at hatch'; bibkey.ah = {'Vigl2023'};
  temp.ah = C2K(15); units.temp.ah = 'K'; label.temp.ah = 'temperature';
  
data.Lh = 0.006852623; units.Lh = 'cm'; label.Lh = 'length at hatch'; bibkey.Lh = {'Vigl2023'};
data.Lb = 0.007812255; units.Lb = 'cm'; label.Lb = 'length at birth'; bibkey.Lb = {'Vigl2023'};
data.Lp = 0.0266569358756218; units.Lp = 'cm'; label.Lp = 'length at puberty'; bibkey.Lp = {'Vigl2023'};

data.Wdp = 0.0000205716902; units.Wdp = 'g'; label.Wdp = 'ultimate dry weight at puberty'; bibkey.Wdp = {'Vigl2023'};

%% set uni-variate data
% time - length
% code modified to work with means and standard deviation
data_tL = readmatrix('tL12.txt');
[t, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL12 = [t, Lmean]; % d, cm- time, mean length
units.tL12 = {'d', 'cm'}; label.tL12 = {'time', 'area^{(1/2)}'};
temp.tL12 = C2K(12); units.temp.tL12 = 'K'; label.temp.tL12 = 'temperature';
bibkey.tL12 = {'Vigl2023'};
stdev.tL12 =Lsd ; units.stdev.tL12 = 'cm'; label.stdev.tL12 = 'standard deviation';
%
data_tL = readmatrix('tL15.txt');
[t, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL15 = [t, Lmean]; % d, cm- time, mean length
units.tL15 = {'d', 'cm'}; label.tL15 = {'time', 'area^{(1/2)}'};
temp.tL15 = C2K(15); units.temp.tL15 = 'K'; label.temp.tL15 = 'temperature';
bibkey.tL15 = {'Vigl2023'};
stdev.tL15 =Lsd ; units.stdev.tL15 = 'cm'; label.stdev.tL15 = 'standard deviation';
%
data_tL = readmatrix('tL18.txt');
[t, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL18 = [t, Lmean]; % d, cm- time, mean length
units.tL18 = {'d', 'cm'}; label.tL18 = {'time', 'area^{(1/2)}'};
temp.tL18 = C2K(18); units.temp.tL18 = 'K'; label.temp.tL18 = 'temperature';
bibkey.tL18 = {'Vigl2023'};
stdev.tL18 =Lsd ; units.stdev.tL18 = 'cm'; label.stdev.tL18 = 'standard deviation';
%
data_tL = readmatrix('tL21.txt');
[t, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL21 = [t, Lmean]; % d, cm- time, mean length
units.tL21 = {'d', 'cm'}; label.tL21 = {'time', 'area^{(1/2)}'};
temp.tL21 = C2K(21); units.temp.tL21 = 'K'; label.temp.tL21 = 'temperature';
bibkey.tL21 = {'Vigl2023'};
stdev.tL21 =Lsd ; units.stdev.tL21 = 'cm'; label.stdev.tL21 = 'standard deviation';
%
data_tL = readmatrix('tL24.txt');
[t, ai, ci] = unique(data_tL (:,1)); maxci = max(ci); 
Lmean = zeros(size(t,1),1); % preallocate zeros for mean L values
Lsd = zeros(size(t,1),1); % preallocate zeros for sd values
for i = 1:maxci
 Lmean(i) =    mean(data_tL (ci==i,2));
 Lsd(i) =    std(data_tL (ci==i,2));
end
data.tL24 = [t, Lmean]; % d, cm- time, mean length
units.tL24 = {'d', 'cm'}; label.tL24 = {'time', 'area^{(1/2)}'};
temp.tL24 = C2K(24); units.temp.tL24 = 'K'; label.temp.tL24 = 'temperature';
bibkey.tL24 = {'Vigl2023'};
stdev.tL24 =Lsd ; units.stdev.tL24 = 'cm'; label.stdev.tL24 = 'standard deviation';

% T-tj data
data.Ttj = [ ... % temperature (C), time since hatch at metam (d)
   14.86   15
   19.34   10
   24.90   7];
  units.Ttj = {'C', 'd'};  label.Ttj = {'temperature', 'time since hatch at metam'};  
  bibkey.Ttj = 'KochBui2017';
  
%%time- number of offspring
data.tN12 = [...
37.00	25
37.00	29
37.00	22
37.00	26
37.00	27
37.00	33
37.00	33
37.00	27
37.00	25
37.00	30
37.00	29
37.00	20];
units.tN12 = {'d', '#'}; label.tN12= {'time', 'clutch size'};
temp.tN12 = C2K(12); units.temp.tN12 = 'K'; label.temp.tN12 = 'temperature';
bibkey.tN12 = {'Vigl2023'}; treat.tN12 = {0};
%
data.tN15 = [...
32.00	41
32.00	25
32.00	10
33.00	31
37.00	33
37.00	31
37.00	21
37.00	25];
units.tN15 = {'d', '#'}; label.tN15= {'time', 'clutch size'};
temp.tN15 = C2K(15); units.temp.tN15 = 'K'; label.temp.tN15 = 'temperature';
bibkey.tN15 = {'Vigl2023'};
treat.tN15 = {0};
%
data.tN18 = [...
23.00	18
24.00	20
24.00	18
24.00	29
24.00	21
30.00	24
30.00	34
30.00	34
30.00	32
30.00	46];
units.tN18 = {'d', '#'}; label.tN18= {'time', 'clutch size'};
temp.tN18 = C2K(18); units.temp.tN18 = 'C'; label.temp.tN18 = 'temperature';
bibkey.tN18 = {'Vigl2023'};
treat.tN18 = {0};
%
data.tN21 = [...
17.00	29
18.00	18
20.00	5
23.00	18
22.00	17
22.00	12
23.00	24
23.00	34
24.00	10];
units.tN21 = {'d', '#'}; label.tN21= {'time', 'clutch size'};
temp.tN21 = C2K(21); units.temp.tN21 = 'K'; label.temp.tN21 = 'temperature';
bibkey.tN21 = {'Vigl2023'};
treat.tN21 = {0};
%
data.tN24 = [...
17.00	12
18.00	21
19.00	10
19.00	23
18.00	19
19.00	10
18.00	11
18.00	12
18.00	11
18.00	10
18.00	16
18.00	15];
units.tN24 = {'d', '#'}; label.tN24= {'time', 'clutch size'};
temp.tN24 = C2K(24); units.temp.tN24 = 'C'; label.temp.tN24 = 'temperature';
bibkey.tN24 = {'Vigl2023'};
treat.tN24 = {0}; 

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


% temperature-max reprod rate
data.TR= [...
12	0.971628449
12	1.127089001
12	0.855033035
12	1.010493587
12	1.049358725
12	1.282549553
12	1.282549553
12	1.049358725
12	0.971628449
12	1.165954139
12	1.127089001
12	0.777302759
15	1.630867144
15	0.994431185
15	0.397772474
15	1.23309467
15	1.312649165
15	1.23309467
15	0.835322196
15	0.994431185
18	0.885391048
18	0.983767831
18	0.885391048
18	1.426463355
18	1.032956222
18	1.180521397
18	1.672405312
18	1.672405312
18	1.574028529
18	2.262666011
21	1.82160804
21	1.130653266
21	0.314070352
21	1.130653266
21	1.067839196
21	0.753768844
21	1.507537688
21	2.135678392
21	0.628140704
24	0.857142857
24	1.5
24	0.714285714
24	1.642857143
24	1.357142857
24	0.714285714
24	0.785714286
24	0.857142857
24	0.785714286
24	0.714285714
24	1.142857143
24	1.071428571];
units.TR = {'C', '#'}; label.TR= {'temperature', 'max reproduction rate'};
bibkey.TR = {'Vigl2023'};

% length-C mass-N mass 
data.LWCN= [...
0.026656936	0.983820513	0.247615385
0.026656937	1.035122449	0.254292517
0.026656938	0.93593141	0.20389359];
% data.LWC = LWCN(:,[1 2]);
units.LWCN = {'cm', 'mugC', 'mugN'}; label.LWCN= {'length', 'carbon mass', 'nitrogen mass'};
bibkey.LWCN = {'Vigl2023'};  
treat.LWCN = {1, {'Carbon weight','Nitrogen weight'}};
% T-tp data
data.Ttp = [ ... % temperature (C), time since metam at puberty (d)
   14.86   12.31
   19.34   10.12
   24.90    8.86];
  units.Ttp = {'C', 'd'};  label.Ttp = {'temperature', 'time since metam at puberty'};  
  bibkey.Ttp = 'KochBui2017';
% temperature age at puberty
data.Tap = [... time since hatch, age at puberty
12  25.73
15  25.14
18  20.33 
21  15.92
24  14];
units.Tap = {'C', 'd'}; label.Tap = {'temperature', 'age at puberty'};
bibkey.Tap = {'Vigl2023'};  

%% set weights for all real data
weights = setweights(data, []);
weights.LWCN= 0* weights.LWCN;

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
auxData.stdev = stdev;
auxData.treat = treat;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
% txtData.comment = comment;

% set1 = {'tL24','tL21', 'tL18','tL15', 'tL12'}; 
% subtitle1 = {'Growth at 24, 21, 18, 15 and 12 C'};
% metaData.grp.sets = {set1, set2};
% metaData.grp.subtitle = {subtitle1,subtitle2};

set1 = {'tN24','tN21', 'tN18','tN15', 'tN12'}; 
subtitle1 = {'Reproduction at 24, 21, 18, 15 and 12 C'};
metaData.grp.sets = {set1};
metaData.grp.subtitle = {subtitle1};


%% Discussion points
D1  = 'There is a change of shape at metamorphosis, between naupliar stage 6 and copepodite 1 stages, and the length is the square root of the ellipse surface';
D2  = 'sbp model is coded as subfunction, and assumes no kappa rule in the adult stage. more investigations pending to resolve the final choice of model for this species.';
metaData.discussion = struct('D1',D1, 'D2',D2);

%% Facts
F1  = 'Sexual reproduction in adult stage; 12 molts: 6 naupliar stages, 5 copepodite stages, 1 adult stage ';
metaData.bibkey.F1 = {'Arnd2013'};
metaData.facts = struct('F1',F1);

%% Links
% metaData.links.id_CoL = '7BDXW'; % Cat of Life
% metaData.links.id_ITIS = '86559'; % ITIS
% metaData.links.id_EoL = '1020216'; % Ency of Life
% metaData.links.id_Wiki = 'Tachidiidae'; % Wikipedia
% metaData.links.id_ADW = 'Tachidius_discipes'; % Anim Div. Web
% metaData.links.id_Taxo = '603913'; % Taxonomicon
% metaData.links.id_WoRMS = '116814'; % WoRMS

%% References
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman, S.A.L.M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press, Cambridge}, ' ...
'pages = {Table 4.2 (page 150), 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Vigl2023'; type = 'misc'; bib = [ ...
'author = {Vigl2023}, ' ... 
'note = {experimental data}, ' ... 
'year = {2023}, ' ... 
'doi = {10.21203/rs.3.rs-2858869/v1}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ', bib, '}'';'];
%
bibkey = 'Arnd2013'; type = 'Phdthesis'; bib = [ ...  
'title = {Testing the suitability of harpacticoid copepods as food for marine fish larvae}, ' ...
'school = {Christian-Albrechts-Universit\"{a}t, Kiel}, ' ...
'year = {2013}, ' ...
'author = {Carmen Arndt}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

