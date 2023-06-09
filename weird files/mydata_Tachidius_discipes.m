function [data, auxData, metaData, txtData, weights] = mydata_Microarthridion_littorale
% file generated by prt_mydata

%% set metaData
metaData.phylum     = 'Arthropoda';
metaData.class      = 'Copepoda';
metaData.order      = 'Harpacticoida';
metaData.family     = 'Tachidiidae';
metaData.species    = 'Microarthridion_littorale';
metaData.species_en = 'no_english_name';

metaData.ecoCode.climate = {'Cfa','Cfb'};
metaData.ecoCode.ecozone = {'MAN','MAn'};
metaData.ecoCode.habitat = {'jiSm','jiM','jiMm','jiMi','jiMcb'};
metaData.ecoCode.embryo  = {'Mbf'};
metaData.ecoCode.migrate = {'TW'};
metaData.ecoCode.food    = {'jiD','jiB','jiPp'};
metaData.ecoCode.gender  = {'D'};
metaData.ecoCode.reprod  = {'O'};

metaData.T_typical  = C2K(15); % K, body temp

metaData.data_0     = {'ah','Lh','Lb','Lp','Li','Wdi','am'};
metaData.data_1     = {'tL','tL','tL','tN'};

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

data.Lh = 0.0082654; units.Lh = 'cm'; label.Lh = 'length at hatch'; bibkey.Lh = {'Vigliano'};
data.Lb = 0.009024; units.Lb = 'cm'; label.Lb = 'length at birth'; bibkey.Lb = {'Vigliano'};
data.Lp = 0.052862; units.Lp = 'cm'; label.Lp = 'length at puberty'; bibkey.Lp = {'Vigliano'};
data.Li = 0.052862; units.Li = 'cm'; label.Li = 'ultimate length'; bibkey.Li = {'Vigliano'};

data.Wdi = 0.0052862; units.Wdi = 'g'; label.Wdi = 'ultimate dry weight'; bibkey.Wdi = {'Vigliano'};

%% set uni-variate data
% time - length
data.tL = [ ... 
  1 0.00861266
  1 0.00796874
  1 0.00806129
  1 0.00812652
  1 0.00979341
  1 0.00805279
  1 0.0080537
  1 0.00790843
  1 0.00807726
  1 0.00815626
  1 0.0077764
  1 0.00850619
  1 0.00823434
  1 0.00828896
  1 0.00798269
  1 0.00800997
  1 0.00796454
  1 0.00844893
  1 0.00830182
  1 0.00807221
  1 0.00806084
  1 0.00855553
  1 0.00794007
  1 0.00778835
  1 0.00809077
  1 0.00809348
  1 0.00786557
  1 0.014116
  1 0.00797185
  1 0.00828759
  1 0.00828963
  1 0.00817336
  1 0.0085632
  1 0.00765103
  1 0.0083357
  1 0.00850633
  1 0.00869011
  1 0.00792964
  1 0.00848126
  1 0.00952285
  1 0.00929334
  1 0.00864166
  1 0.0080043
  1 0.00804064
  1 0.00789665
  1 0.0078
  1 0.00808775
  1 0.00821678
  1 0.00814713
  1 0.00827009
  2 0.0090248
  2 0.00834064
  2 0.00790265
  2 0.00892798
  2 0.00788482
  2 0.00770743
  2 0.00802953
  2 0.0084855
  1 0.00837902
  1 0.00800023
  1 0.00797205
  1 0.00790019
  1 0.00836016
  1 0.00787076
  1 0.00865922
  1 0.00753329
  1 0.00764881
  1 0.00859443
  1 0.00891099
  1 0.00778105
  1 0.00814915
  1 0.00824025
  1 0.00876176
  1 0.00848928
  1 0.00770349
  1 0.0088901
  1 0.0079939
  1 0.00817531
  1 0.00817654
  2 0.00810337
  2 0.00846806
  2 0.00853651
  2 0.00898986
  2 0.00853012
  2 0.00882149
  2 0.00848671
  2 0.00993185
  2 0.00897602
  2 0.00875456
  2 0.00897863
  2 0.00936851
  2 0.00757399
  2 0.00921684
  2 0.00952942
  2 0.00893776
  2 0.00829591
  2 0.00889659
  2 0.00897074
  2 0.0109022
  2 0.00916384
  2 0.00816434
  2 0.00993532
  2 0.00817711
  2 0.00859211
  2 0.0104803
  2 0.0107156
  2 0.00846232
  2 0.00794822
  2 0.00779272
  3 0.00847768
  3 0.00884446
  3 0.00844005
  3 0.00970763
  2 0.00951298
  2 0.00942869
  2 0.00931834
  2 0.00929418
  2 0.00838999
  2 0.00956955
  2 0.00781427
  2 0.00905868
  2 0.0094774
  2 0.00896712
  2 0.00965205
  2 0.00868947
  2 0.00975618
  2 0.00828022
  2 0.00980592
  2 0.00943486
  2 0.00923632
  2 0.00985014
  2 0.00998482
  2 0.00959698
  3 0.00999061
  3 0.00909369
  3 0.00914899
  3 0.00965989
  3 0.00997557
  3 0.010094
  3 0.00927396
  3 0.00805283
  3 0.0101829
  3 0.00993295
  3 0.00908103
  3 0.0111552
  3 0.0105705
  4 0.0103236
  4 0.0109543
  4 0.0147868
  4 0.0121813
  4 0.0123623
  4 0.0114791
  4 0.00952065
  4 0.012502
  4 0.0104125
  4 0.0118705
  4 0.0120717
  4 0.00877933
  4 0.00827225
  4 0.0105418
  4 0.0110052
  4 0.00817429
  4 0.0100591
  4 0.01217
  4 0.0102786
  4 0.0134112
  4 0.00883464
  4 0.0108013
  4 0.00919247
  5 0.00982914
  5 0.0115444
  5 0.0112525
  5 0.00921071
  5 0.0116669
  5 0.00918784
  5 0.012739
  5 0.00916249
  5 0.0145703
  5 0.00986672
  4 0.0120139
  4 0.0106145
  4 0.0108408
  4 0.0109089
  4 0.0104601
  4 0.0110738
  4 0.010941
  4 0.010943
  4 0.00963355
  4 0.0100173
  4 0.0133938
  4 0.00922599
  4 0.0132887
  4 0.0121317
  4 0.0101932
  4 0.0114294
  4 0.0116111
  4 0.0130359
  4 0.0110162
  4 0.0117742
  4 0.011575
  5 0.0118434
  5 0.0115445
  5 0.0109051
  5 0.0129142
  5 0.0113646
  5 0.0114382
  5 0.011514
  5 0.0109092
  5 0.0112439
  5 0.0108916
  5 0.0134883
  5 0.0138362
  5 0.0137077
  5 0.0148471
  5 0.0111344
  5 0.0102163
  5 0.0112173
  5 0.0106821
  5 0.0108281
  5 0.0115282
  5 0.0113466
  5 0.0130025
  5 0.0134912
  5 0.0116654
  5 0.0124334
  5 0.0140888
  5 0.0109468
  5 0.013536
  5 0.0124168
  5 0.0152921
  5 0.012614
  5 0.0122664
  5 0.0126486
  5 0.0114302
  5 0.0121528
  5 0.0158376
  5 0.0138664
  5 0.0117217
  5 0.0109138
  5 0.0131065
  5 0.0124454
  5 0.0148009
  5 0.01295
  6 0.0160347
  6 0.010952
  6 0.0150957
  6 0.0106771
  6 0.0158785
  5 0.0131385
  5 0.0136918
  5 0.0157603
  5 0.0103915
  5 0.0147028
  5 0.0117846
  5 0.0112156
  5 0.0141056
  5 0.0130495
  5 0.0162445
  6 0.0154205
  6 0.0158351
  6 0.0183292
  6 0.0166235
  6 0.0154346
  6 0.0149927
  6 0.0142185
  6 0.013593
  6 0.0135953
  6 0.0163568
  6 0.0140165
  6 0.0136222
  6 0.0144381
  6 0.0140291
  6 0.0147432
  6 0.0153998
  6 0.015387
  6 0.013926
  6 0.0149602
  6 0.0151043
  6 0.0148981
  6 0.0142153
  6 0.0123396
  6 0.0135232
  6 0.0146713
  6 0.014138
  6 0.0141219
  6 0.0142572
  6 0.0146353
  6 0.0108721
  6 0.015695
  6 0.0144359
  6 0.0157769
  6 0.0144024
  6 0.0126775
  6 0.0164239
  6 0.0139856
  6 0.0133351
  6 0.0148174
  6 0.0151206
  6 0.0120581
  6 0.0128273
  6 0.0161308
  6 0.0121538
  6 0.0141332
  6 0.0157032
  6 0.022009
  6 0.0151577
  6 0.0138639
  6 0.0121135
  7 0.00946285
  7 0.0104343
  7 0.0188658
  7 0.0241694
  7 0.010194
  7 0.0121392
  7 0.0166393
  6 0.0141554
  6 0.0151776
  6 0.0142468
  6 0.0149621
  6 0.0152118
  6 0.0160967
  6 0.0165857
  6 0.0160759
  6 0.010981
  6 0.016801
  6 0.0152618
  6 0.0180562
  6 0.00955548
  6 0.0168998
  6 0.0157756
  6 0.0175788
  6 0.0120623
  6 0.0152523
  6 0.0167622
  7 0.0149111
  7 0.0159753
  7 0.0203931
  7 0.0178151
  7 0.0176736
  7 0.0148801
  7 0.0166631
  7 0.0166339
  7 0.0177414
  7 0.018337
  7 0.0184
  7 0.0179794
  7 0.0182336
  7 0.0180504
  7 0.0180802
  7 0.016892
  7 0.0182873
  7 0.0188747
  7 0.0238062
  7 0.0170141
  7 0.0174109
  7 0.0161752
  7 0.0176389
  7 0.0165934
  7 0.0168859
  7 0.0119699
  7 0.0268686
  7 0.0159092
  7 0.0166036
  7 0.0188752
  7 0.0141093
  7 0.0146112
  7 0.0145205
  8 0.0157761
  8 0.0168185
  8 0.0125709
  8 0.0105069
  8 0.0184711
  8 0.0175521
  7 0.0155712
  7 0.0211406
  7 0.0157615
  7 0.0178243
  7 0.0167913
  7 0.0156379
  7 0.0159179
  7 0.015571
  7 0.0165326
  7 0.0172324
  7 0.0198296
  7 0.0114291
  7 0.0173924
  7 0.0167976
  7 0.0168599
  7 0.0157761
  8 0.0174559
  8 0.0150312
  8 0.0153684
  8 0.0138654
  8 0.0161042
  8 0.0157285
  8 0.0155231
  8 0.0248131
  8 0.0175726
  8 0.016907
  8 0.017843
  8 0.014671
  8 0.0142361
  8 0.017455
  8 0.0171346
  8 0.0166708
  8 0.0138908
  8 0.0141616
  8 0.0158978
  8 0.0171334
  8 0.0180172
  8 0.0149048
  8 0.0159369
  8 0.0161656
  8 0.0133614
  8 0.0125053
  8 0.0161425
  8 0.0141325
  8 0.0154462
  8 0.0173205
  8 0.0164241
  8 0.0155837
  8 0.0163216
  9 0.0163171
  9 0.0155831
  9 0.016897
  9 0.018168
  8 0.0175066
  8 0.015344
  8 0.0167604
  8 0.017111
  8 0.0168216
  8 0.0165183
  8 0.0172003
  8 0.0179768
  8 0.0166001
  8 0.0166303
  8 0.0149202
  8 0.0167526
  8 0.0101133
  8 0.0167252
  8 0.0169115
  8 0.0147692
  9 0.0154104
  9 0.0178195
  9 0.016388
  9 0.0178472
  9 0.0167807
  9 0.01632
  9 0.0178699
  9 0.0226482
  9 0.0154362
  9 0.019129
  9 0.0182299
  9 0.0146333
  9 0.0175157
  9 0.0155046
  9 0.0217765
  10 0.0263795
  10 0.0261932
  10 0.0266951
  10 0.02918
  10 0.0251029
  10 0.0149291
  10 0.0182307
  10 0.0254348
  10 0.02107
  10 0.0208691
  10 0.0240182
  10 0.0201641
  10 0.018363
  10 0.0263164
  10 0.0262177
  10 0.0249831
  10 0.0244578
  10 0.0257634
  11 0.0251204
  11 0.0300422
  11 0.0268781
  11 0.0198709
  11 0.0301481
  10 0.0213647
  10 0.0199257
  10 0.0246605
  10 0.0253517
  10 0.0252161
  10 0.0198721
  10 0.0258699
  10 0.0262479
  10 0.0206181
  10 0.0249528
  10 0.0240642
  10 0.0240014
  10 0.0125087
  10 0.0224115
  10 0.0248297
  10 0.0249957
  11 0.0234649
  11 0.0256819
  11 0.0236338
  11 0.0242652
  11 0.0266087
  11 0.0269291
  11 0.0264679
  11 0.0240343
  11 0.0307172
  11 0.0311535
  11 0.030657
  11 0.0312378
  11 0.0321584
  11 0.0304034
  11 0.0280668
  11 0.0253129
  11 0.0310184
  11 0.0333164
  11 0.0364098
  11 0.0253963
  11 0.0243462
  11 0.0287959
  11 0.0290179
  11 0.0228505
  11 0.0241517
  11 0.0238395
  11 0.03133
  12 0.02522
  12 0.0304588
  12 0.0243667
  12 0.0307696
  11 0.0254395
  11 0.02465
  11 0.0239935
  11 0.0277627
  11 0.0251315
  11 0.0241623
  11 0.024813
  11 0.0246604
  11 0.0214599
  11 0.0244259
  11 0.0230199
  11 0.0249456
  12 0.0259041
  12 0.0296699
  12 0.0218472
  12 0.0239003
  12 0.0279921
  12 0.0254479
  12 0.0324559
  12 0.03468
  12 0.0314141
  12 0.0300131
  12 0.0310212
  12 0.0110727
  12 0.0287681
  12 0.0292101
  12 0.0268204
  12 0.0163096
  12 0.0307165
  12 0.0267854
  12 0.0253463
  12 0.0362939
  12 0.0135749
  12 0.0318077
  12 0.0297931
  12 0.0244005
  12 0.0247129
  13 0.0311826
  13 0.0245063
  13 0.020691
  13 0.0340763
  13 0.0280885
  12 0.0273391
  12 0.0295791
  12 0.0211819
  12 0.0239727
  12 0.0296575
  12 0.0233454
  12 0.0237233
  12 0.0215457
  12 0.0242631
  12 0.0244059
  13 0.0237652
  13 0.0264451
  13 0.0276269
  13 0.0369635
  13 0.0326444
  13 0.0305127
  13 0.0355029
  13 0.0302231
  13 0.0355721
  13 0.0367506
  13 0.026453
  13 0.0224485
  13 0.0276735
  13 0.0263323
  13 0.0259004
  13 0.0243082
  14 0.031968
  14 0.0279959
  14 0.0354925
  13 0.0237123
  13 0.0259784
  13 0.0228113
  13 0.0283672
  13 0.0285657
  13 0.0250905
  13 0.0264815
  13 0.0209538
  13 0.0246449
  13 0.0118445
  13 0.0202982
  13 0.0213342
  13 0.0245816
  13 0.0237822
  13 0.0255437
  13 0.020878
  14 0.0286155
  14 0.0236116
  14 0.0287186
  14 0.0327521
  14 0.0260139
  14 0.0276141
  14 0.0366709
  14 0.0379572
  14 0.0329324
  14 0.0376851
  14 0.0364423
  14 0.0390535
  14 0.0355037
  15 0.0330092
  15 0.0323915
  15 0.0335113
  15 0.0173647
  15 0.0336378
  15 0.0336031
  15 0.0293789
  14 0.0322955
  14 0.0335393
  14 0.0293709
  14 0.0309608
  14 0.0299798
  14 0.0263335
  14 0.0335975
  14 0.0345961
  15 0.0434288
  15 0.0359669
  15 0.0423371
  15 0.0358376
  15 0.04559
  15 0.0348838
  15 0.0447662
  14 0.0341996
  14 0.036911
  14 0.0277127
  14 0.026244
  14 0.0283859
  14 0.0283556
  14 0.0298851
  14 0.0281408
  14 0.0204862
  14 0.0231238
  14 0.0216585
  14 0.0226657
  14 0.0224907
  14 0.0235984
  14 0.025187
  15 0.0398934
  15 0.0352969
  15 0.0384689
  15 0.039226
  15 0.0428481
  15 0.0385925
  15 0.0428745
  15 0.0396686
  15 0.0365379
  15 0.0356861
  15 0.0330082
  15 0.033588
  15 0.0409104
  15 0.0381992
  15 0.0349633
  15 0.0377033
  15 0.0338326
  15 0.0328128
  15 0.0388333
  15 0.0275952
  15 0.037255
  15 0.0391471
  15 0.0380049
  15 0.0335479
  16 0.0328938
  16 0.0272885
  16 0.035755
  16 0.0356509
  15 0.024
  15 0.0288418
  15 0.0247468
  15 0.0191157
  15 0.0292947
  15 0.0200872
  15 0.0221594
  15 0.0229878
  16 0.0257593
  16 0.0276993
  16 0.0273353
  16 0.0214205
  16 0.0262109
  16 0.0245176
  16 0.0399191
  16 0.0402409
  16 0.0410167
  16 0.042301
  16 0.0436915
  16 0.0100986
  16 0.0357145
  16 0.0369094
  16 0.035834
  16 0.0330784
  16 0.0410345
  16 0.0373006
  16 0.0317814
  16 0.0347077
  16 0.0362865
  16 0.0349689
  17 0.0555346
  17 0.0336781
  17 0.0294015
  17 0.0422448
  17 0.0214144
  17 0.0407151
  17 0.0364172
  17 0.0389765
  17 0.0378204
  17 0.0321555
  16 0.0251428
  16 0.0245229
  16 0.025086
  16 0.0237082
  16 0.024916
  16 0.0240468
  17 0.0285734
  17 0.0281273
  17 0.022685
  17 0.0264161
  17 0.0303889
  17 0.024061
  17 0.040185
  17 0.0435943
  17 0.041599
  17 0.0340869
  17 0.0410897
  17 0.0408727
  17 0.0412242
  17 0.045976
  18 0.0439843
  18 0.0382371
  18 0.0370096
  18 0.0384998
  18 0.0331408
  18 0.0286353
  18 0.02989
  18 0.0304027
  18 0.0360033
  18 0.0334603
  18 0.0305981
  19 0.0342738
  19 0.042386
  19 0.0487445
  19 0.0349827
  19 0.0427475
  19 0.0432064
  18 0.028762
  18 0.0274296
  18 0.0302864
  18 0.0311452
  18 0.0216632
  18 0.0277658
  18 0.0223445
  18 0.0231896
  18 0.0314463
  18 0.036881
  19 0.0285591
  19 0.0283719
  19 0.0292887
  19 0.0294164
  19 0.0307664
  19 0.0438239
  19 0.0418084
  19 0.0349015
  19 0.0455361
  19 0.0444103
  19 0.0502143
  19 0.034892
  19 0.0374529
  19 0.0487297
  19 0.0371487
  19 0.0409174
  19 0.0348823
  19 0.0360655
  19 0.0279615
  19 0.0322898
  19 0.0330539
  19 0.0302616
  20 0.047704
  20 0.0343726
  20 0.0455337
  20 0.0440621
  20 0.0532663
  19 0.0325989
  19 0.0263702
  19 0.0345492
  19 0.0316279
  19 0.0325222
  19 0.0344551
  19 0.024416
  19 0.0252702
  19 0.0266029
  19 0.0286469
  19 0.0239319
  20 0.0255401
  20 0.0322673
  20 0.0344729
  20 0.0297577
  20 0.0276207
  20 0.0370504
  20 0.0425431
  20 0.0434708
  20 0.0388216
  20 0.0422053
  20 0.045134
  20 0.037778
  20 0.0504505
  20 0.0511437
  20 0.0596848
  20 0.0462703
  20 0.0577474
  20 0.0572134
  20 0.0593791
  20 0.046
  20 0.0553307
  20 0.050414
  20 0.0458718
  20 0.0550326
  21 0.0595039
  21 0.0694158
  21 0.0621999
  21 0.0411696
  21 0.0570912
  20 0.0435867
  20 0.0457824
  20 0.0467632
  20 0.033796
  20 0.0406402
  20 0.028356
  20 0.0333802
  20 0.0303488
  20 0.0280662
  20 0.0337003
  21 0.0336522
  21 0.0418699
  21 0.0351639
  21 0.0464185
  21 0.0457793
  21 0.0567042
  21 0.0577808
  21 0.0548479
  21 0.0506545
  21 0.0616781
  21 0.0649357
  21 0.041995
  21 0.0370391
  21 0.0466662
  21 0.0441074
  21 0.0363389
  21 0.0405104
  21 0.0330117
  21 0.0317274
  21 0.0361792
  21 0.0324085
  21 0.0319606
  22 0.0452437
  22 0.0511127
  22 0.049827
  22 0.0520063
  22 0.0473007
  21 0.0348411
  21 0.0303948
  21 0.0280924
  21 0.0321488
  21 0.0394348
  21 0.0398223
  21 0.0255069
  21 0.0413308
  21 0.0374439
  21 0.0353974
  21 0.0243882
  21 0.0516482
  21 0.0326464
  22 0.0320421
  22 0.0295803
  22 0.0288815
  22 0.0281305
  22 0.0310384
  22 0.0313517
  22 0.0475324
  22 0.0429161
  22 0.0514693
  22 0.0447897
  22 0.0496179
  22 0.0449217
  22 0.0607467
  22 0.0441583
  22 0.0407089
  22 0.0481034
  22 0.0481997
  22 0.0479927
  22 0.0460167
  22 0.0295946
  22 0.0317973
  22 0.032656
  22 0.0275347
  22 0.0291459
  23 0.0473808
  23 0.0445303
  23 0.0485416
  23 0.0503333
  23 0.0554038
  22 0.0356269
  22 0.0303913
  22 0.0355242
  22 0.0291641
  22 0.0335135
  22 0.0327472
  22 0.0339724
  22 0.02576
  22 0.0432718
  22 0.0318173
  22 0.0254306
  23 0.0345327
  23 0.0289104
  23 0.028158
  23 0.035399
  23 0.0347097
  23 0.0408757
  23 0.0466506
  23 0.0388152
  23 0.050452
  23 0.043528
  23 0.0489883
  23 0.0460231
  23 0.0496577
  23 0.0469885
  23 0.0513669
  23 0.0459038
  23 0.0355425
  23 0.0339337
  23 0.0358806
  23 0.0286405
  23 0.0373439
  23 0.0311264
  23 0.0332008
  24 0.0534874
  24 0.040069
  24 0.0384613
  24 0.0585594
  24 0.0536502
  23 0.0354843
  23 0.0309978
  23 0.0368046
  23 0.0397632
  23 0.0399191
  23 0.039432
  23 0.0390311
  23 0.0384209
  23 0.0368451
  23 0.0572082
  23 0.0280043
  23 0.0266502
  23 0.0284599
  23 0.0277474
  24 0.0376453
  24 0.0406264
  24 0.037205
  24 0.0380365
  24 0.0352892
  24 0.0355284
  24 0.0350051
  24 0.0364866
  24 0.0333994
  24 0.0459029
  24 0.0448039
  24 0.0516059
  24 0.0296786
  24 0.0507508
  24 0.0489186
  24 0.0463514
  24 0.043103
  24 0.0417032
  24 0.0527086
  24 0.0481582
  24 0.0595907
  24 0.045632
  24 0.0397845
  24 0.0371239
  24 0.035426
  24 0.0419938
  24 0.0411636
  24 0.0385891
  24 0.0476804
  24 0.0376241
  24 0.0408128
  26 0.0450898
  26 0.0442047
  26 0.0505964
  26 0.0538783
  26 0.0549607
  26 0.0585841
  26 0.0572369
  24 0.0430752
  24 0.0376745
  24 0.0457704
  24 0.0396964
  24 0.038019
  24 0.0441019
  24 0.0358155
  24 0.0431937
  24 0.0442367
  24 0.0416261
  24 0.0272956
  24 0.0314168
  24 0.0331942
  24 0.0324136
  24 0.0363056
  26 0.0349206
  26 0.040677
  26 0.0389842
  26 0.0362087
  26 0.0281645
  26 0.039974
  26 0.0446984
  26 0.0436122
  26 0.0506991
  26 0.0474044
  26 0.0496202
  26 0.0442186
  26 0.0458062
  26 0.0505553
  26 0.0506645
  26 0.0521457
  26 0.0474947
  26 0.0522298
  26 0.0486354
  26 0.0406205
  26 0.0460162
  26 0.0459503
  26 0.0473169
  25 0.0486683
  25 0.0509544
  25 0.0463402
  25 0.0591918
  25 0.0545446
  25 0.0519936
  25 0.0375869
  25 0.0320316
  25 0.0415878
  25 0.0437952
  25 0.044158
  25 0.0432128
  25 0.0432947
  25 0.0488793
  25 0.0475187
  27 0.0486263
  27 0.0548758
  27 0.055106
  25 0.0507281
  25 0.0395172
  25 0.0399948
  25 0.0346308
  25 0.0407078
  25 0.0432797
  25 0.0430302
  25 0.0266467
  25 0.0424007
  25 0.0406107
  25 0.033202
  25 0.0402531
  25 0.0285819
  25 0.0422696
  27 0.038321
  27 0.0398176
  27 0.040254
  27 0.0350827
  27 0.0377703
  27 0.0415629
  27 0.0489625
  27 0.0306543
  27 0.0472957
  27 0.0487526
  27 0.0526052
  27 0.0504585
  27 0.0483887
  27 0.0533496
  27 0.0493828
  27 0.0492704
  27 0.0513087
  26 0.0578173
  26 0.0417612
  26 0.0420251
  26 0.0358645
  26 0.0447106
  26 0.0499486
  26 0.0478261
  26 0.0445889
  26 0.0450165
  26 0.044232
  26 0.0472395
  28 0.0523282
  28 0.0598791
  26 0.0363087
  26 0.0464091
  26 0.0421043
  26 0.0459609
  26 0.0460052
  26 0.0482826
  26 0.0512437
  26 0.0476144
  26 0.0396982
  26 0.0371175
  26 0.0317874
  26 0.0453798
  26 0.0374711
  28 0.0391085
  28 0.0434747
  28 0.0394752
  28 0.0463198
  28 0.0420621
  28 0.0395336
  27 0.0455855
  27 0.0586139
  27 0.0437088
  27 0.0565237
  27 0.0537895
  27 0.0535278
  27 0.0413913
  27 0.0384776
  27 0.0429443
  27 0.0395051
  27 0.0363489
  29 0.0442167
  29 0.0497433
  29 0.0500892
  29 0.0511378
  29 0.0510174
  29 0.0555886
  29 0.0445002
  29 0.0390837
  29 0.0396922];
units.tL = {'d', 'cm'}; label.tL = {'time', 'length'};
temp.tL = C2K(15); units.temp.tL = 'K'; label.temp.tL = 'temperature';
bibkey.tL = {'Vigliano'};
comment.tL = 'experimental data';

%% set weights for all real data
weights = setweights(data, []);

%% set pseudodata and respective weights
[data, units, label, weights] = addpseudodata(data, units, label, weights);

%% pack auxData and txtData for output
auxData.temp = temp;
txtData.units = units;
txtData.label = label;
txtData.bibkey = bibkey;
txtData.comment = comment;

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
