function [data, auxData, metaData, txtData, weights] = mydata_Tachidius_discipes
% file generated by prt_mydata

%% set metaData
metaData.phylum     = 'Arthropoda';
metaData.class      = 'Copepoda';
metaData.order      = 'Harpacticoida';
metaData.family     = 'Tachidiidae';
metaData.species    = 'Tachidius_discipes';
metaData.species_en = 'no_english_name';

metaData.ecoCode.climate = {'Dfb','Cfb'};
metaData.ecoCode.ecozone = {'MAb','MAn'};
metaData.ecoCode.habitat = {'jiSm','jiM','jiMm','jiMi','jiMcb'};
metaData.ecoCode.embryo  = {'Mbf'};
metaData.ecoCode.migrate = {'TW'};
metaData.ecoCode.food    = {'jiD','jiB','jiPp'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(15); % K, body temp

metaData.data_0     = {'ah','Lh','Lb','Lp','Li','Wdi','am'};
metaData.data_1     = {'tL','tL','tL'};

metaData.COMPLETE   = 2; % using criteria of LikaKear2011

metaData.author     = {'Julieta Vigliano Relva'};
metaData.date_subm  = [2023 5 22];
metaData.email      = {'julieta.vigliano@ugent.be'};
metaData.address    = {'UGent Ghent University'};

metaData.curator    = {'Starrlight Augustine'};
metaData.email_cur  = {'starrlight@tecnico.ulisboa.pt'};
metaData.date_acc   = [2023 5 3];

%% set zero-variate data
data.ah = 4.375; units.ah = 'd'; label.ah = 'age at hatch'; bibkey.ah = {'Vigliano'};
  temp.ah = C2K(15); units.temp.ah = 'K'; label.temp.ah = 'temperature';
data.Lh = 0.006852623; units.Lh = 'cm'; label.Lh = 'length at hatch'; bibkey.Lh = {'Vigliano'};
data.Lb = 0.007812255; units.Lb = 'cm'; label.Lb = 'length at birth'; bibkey.Lb = {'Vigliano'};
data.Lp = 0.0266569358756218; units.Lp = 'cm'; label.Lp = 'length at puberty'; bibkey.Lp = {'Vigliano'};
data.Wdp =0.0000205716902; units.Wdp = 'g'; label.Wdp = 'ultimate dry weight at puberty'; bibkey.Wdp = {'Vigliano'};

%% set uni-variate data
% time - length
data.tL12 = readmatrix('tL12.txt');
units.tL12 = {'d', 'cm'}; label.tL12 = {'time', 'length'};
temp.tL12 = C2K(12); units.temp.tL12 = 'K'; label.temp.tL12 = 'temperature';
bibkey.tL12 = {'Vigliano'};

data.tL15 = readmatrix('tL15.txt');
units.tL15 = {'d', 'cm'}; label.tL15 = {'time', 'length'};
temp.tL15 = C2K(15); units.temp.tL = 'K'; label.temp.tL15 = 'temperature';
bibkey.tL15 = {'Vigliano'};

data.tL18 = readmatrix('tL18.txt');
units.tL18 = {'d', 'cm'}; label.tL18 = {'time', 'length'};
temp.tL18 = C2K(18); units.temp.tL18 = 'K'; label.temp.tL18 = 'temperature';
bibkey.tL18 = {'Vigliano'};

data.tL21= readmatrix('tL21.txt');
units.tL21 = {'d', 'cm'}; label.tL21 = {'time', 'length'};
temp.tL21 = C2K(21); units.temp.tL21 = 'K'; label.temp.tL21 = 'temperature';
bibkey.tL21 = {'Vigliano'};

data.tL24 = readmatrix('tL21.txt');
units.tL24 = {'d', 'cm'}; label.tL24 = {'time', 'length'};
temp.tL24 = C2K(24); units.temp.tL24 = 'K'; label.temp.tL24 = 'temperature';
bibkey.tL24 = {'Vigliano'};

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
37.00	20
];
units.tN12 = {'d', '#'}; label.tN12= {'time', 'clutch size'};
temp.tN12 = C2K(12); units.temp.tN12 = 'K'; label.temp.tN12 = 'temperature';
bibkey.tL24 = {'Vigliano'}; treat.tN12 = {0};


data.tN15 = [...
32.00	41
32.00	25
32.00	10
33.00	31
37.00	33
37.00	31
37.00	21
37.00	25
];
units.tN15 = {'d', '#'}; label.tN15= {'time', 'clutch size'};
temp.tN15 = C2K(15); units.temp.tN15 = 'K'; label.temp.tN15 = 'temperature';
bibkey.tL15 = {'Vigliano'};
treat.tN15 = {0};

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
30.00	46
];
units.tN18 = {'d', '#'}; label.tN18= {'time', 'clutch size'};
temp.tN18 = C2K(18); units.temp.tN18 = 'K'; label.temp.tN18 = 'temperature';
bibkey.tN18 = {'Vigliano'};
treat.tN18 = {0};

data.tN21 = [...
17.00	29
18.00	18
20.00	5
23.00	18
22.00	17
22.00	12
23.00	24
23.00	34
24.00	10
];
units.tN21 = {'d', '#'}; label.tN21= {'time', 'clutch size'};
temp.tN21 = C2K(21); units.temp.tL21 = 'K'; label.temp.tN21 = 'temperature';
bibkey.tN_21 = {'Vigliano'};
treat.tN21 = {0};

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
18.00	15
];
units.tN24 = {'d', '#'}; label.tN24= {'time', 'clutch size'};
temp.tN24 = C2K(24); units.temp.tN24 = 'K'; label.temp.tN24 = 'temperature';
bibkey.tN24 = {'Vigliano'};
treat.tN24 = {0}; 

%%tN_T data
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
units.Tah = {'K', 'days'}; label.Tah= {'temperature', 'time'};
bibkey.Tah = {'Vigliano'};

data.LWC=[...
0.026656936	1.009471481
0.026656937	0.983820513
0.026656938	1.035122449
];
units.LWC = {'cm', 'mug'}; label.LWC= {'length', 'carbon mass'};
bibkey.LWC = {'Vigliano'};

%% set weights for all real data
weights = setweights(data, []);
weights.LWC = 10* weights.LWC; 

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
auxData.treat = treat;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
% txtData.comment = comment;

set1 = {'tL12','tL15', 'tL18','tL21', 'tL24'}; 
subtitle1 = {'Growth at different temperatures'};
set2 = {'tN12','tN15', 'tN18','tN21', 'tN24'}; 
subtitle2 = {'Reproduction at 12, 15, 18, 21 and 24 debC'};
metaData.grp.sets = {set1, set2};
metaData.grp.comment = {subtitle1,subtitle2};

%% Discussion points
D1  = 'There is a change of shape at metamorphosis, between naupliar stage 6 and copepodite 1 stages';
% metaData.bibkey.D1 = {};
metaData.discussion = struct('D1',D1);

%% Facts
F1  = 'Sexual reproduction in adult stage; 12 molts: 6 naupliar stages, 5 copepodite stages, 1 adult stage ';
% metaData.bibkey.F1 = {};
metaData.facts = struct('F1',F1);

%% Links
metaData.links.id_CoL = '42RK9'; % Cat of Life
metaData.links.id_EoL = '46538747'; % Ency of Life
metaData.links.id_ADW = 'Microarthridion_littorale'; % Anim Div. Web
metaData.links.id_Taxo = '465425'; % Taxonomicon
metaData.links.id_WoRMS = '116498';


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

bibkey = 'Vigliano'; type = 'misc'; bib = [ ...
'author = {Vigliano}, ' ... 
'note = {experimental data}, ' ... 
'year = {2023}, ' ... 
'doi = {10.21203/rs.3.rs-2858869/v1}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ', bib, '}'';'];
%


