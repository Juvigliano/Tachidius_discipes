%% fig:altri_gallus
%% bib:RomiLokh51
%% out:romi51w romi51o
%% Gallus domesticus: wet weight and dioxygen consumption at age

aW_Gd = [4.922633708 0.2214927143;
      5.902508118 0.6859380295;
      6.941072921 0.6301200864;
      6.941886525 0.2061112758;
      6.959415981 1.092712216;
      7.881857383 1.422134446;
      7.883188734 0.7283018468;
      7.884039319 0.2850199085;
      7.900422334 1.769087809;
      9.879142906 2.832927397;
      9.918602673 2.312627076;
      9.919712132 1.734433243;
      10.87680564 4.049061759;
      10.87761925 3.625052948;
      11.97317328 3.935535624;
      13.0457615  6.192567138;
      13.04757361 5.248183878;
      13.17638183 8.274324382;
      14.09582769 10.16486996;
      14.11147107 12.03440042;
      15.02855008 15.15842617;
      15.03642725 11.05324996;
      15.90572555 19.03402978;
      15.99011842 15.14100623;
      15.99525891 12.46204147;
      16.88452748 20.0573958;
      16.94037026 21.02116338;
      17.99210065 24.12617545;
      18.0334465  22.62294561;
      19.96046627 30.58649081];

aO2_Gd = [2.106267042 5.114567709;
       3.141336813 5.032435329;
       3.931421041 8.138761959;
       5.089711856 10.15953171;
       5.935270479 27.69810086;
       6.936779043 35.36512367;
       7.921507003 46.90672407;
       8.977172137 62.31594169;
       9.947963222 92.16853843;
       10.99496278 152.6489432;
       11.95986267 241.6570382;
       12.82109648 338.7719643;
       13.84618719 466.1556304;
       14.79735394 485.446386;
       15.84289696 538.5325276;
       16.84822028 565.565478;
       17.86470611 560.203121;
       18.86108209 541.8141685]; 

%% we must have decreasing times for embryo-data fitting
aW_Gd = aW_Gd(30:-1:1,:); aO2_Gd = aO2_Gd(18:-1:1,:);
%% avoid large numbers for numerical stability
aO2_Gd(:,2) = .001 * aO2_Gd(:,2); 
aO2_Gd = [aO2_Gd, 10 * aO2_Gd(:,2)]; aO2_Gd(11:18,3) = 0;
