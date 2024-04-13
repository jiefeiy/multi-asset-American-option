clear, clc, close all
% reference d = 2:16
load('v_ref.mat')

% using SGGK with qlevel = 4;
dim = 2:12;
ilevel = 3:7;

nb_dim = size(dim,2);
nb_ilevel = size(ilevel,2);

value = cell(1,nb_dim);
nb_knots_all = cell(1,nb_dim);

% d = 2
value{1} = [2.948949116402245   3.187983360746620   3.188038515370521   3.183899095187627   3.183145582488622];
nb_knots_all{1} = [29    13   173;
    65    33   173;
   145    81   173;
   321   193   173;
   705   449   173];

% d = 3
value{2} = [2.912972628181152   3.022587278584994   3.005981682951669   3.002939802264737   3.003015785625988];
nb_knots_all{2} = [69    19   495;
    177    55   495;
    441   151   495;
    1073   399   495;
    2561   1023   495];

% d = 4
value{3} = [2.902798711967551   2.931449889059475   2.910323821975349   2.908773749642587];
nb_knots_all{3} = [137          25        1097;
         401          81        1097;
        1105         241        1097;
        2929         689        1097];

% d = 5
value{4} = [2.896267142532981   2.871594363928806   2.851374056702951   2.850511276329679];
nb_knots_all{4} = [241          31        2091;
    801         111        2091;
    2433         351        2091;
    6993        1071        2091];

% d = 6
value{5} = [2.888908545655187   2.828251980507589   2.811030131073752   2.810723588593970];
nb_knots_all{5} = [389          37        3605;
    1457         145        3605;
    4865         481        3605;
    15121        1553        3605];

% d = 7
value{6} = [2.880974824757323   2.795528612778188   2.781687635869414   2.781742763955954];
nb_knots_all{6} = [589          43        5783;
    2465         183        5783
    9017         631        5783;
    30241        2143        5783];

% d = 8
value{7} = [2.872484022528131   2.771787292212994   2.759895636709271   2.760205084789753];
nb_knots_all{7} = [849          49        8785;
    3937         225        8785;
    15713         801        8785;
    56737        2849        8785];

% d = 9
value{8} = [2.863102002121884   2.750536030751086   2.742842295394100];
nb_knots_all{8} = [1177          55       12787;
    6001         271       12787;
    26017         991       12787];

% d = 10
value{9} = [2.852525075164261   2.731982479479699   2.729221845867784];
nb_knots_all{9} = [1581	  61	  17981
8801	321	   17981
41265	1201	17981];

% d = 11
value{10} = [2.841182963097624   2.715370997105701   2.7180];
nb_knots_all{10} = [2069          67       24575;
       12497         375       24575;
       63097        1431       24575];

% d = 12
value{11} = [2.829357375299730   2.699951745447553   2.708529647502717];
nb_knots_all{11} = [2649          73       32793;
       17265         433       32793;
       93489        1681       32793];

disp('-----------reference price-----------');
disp(vpa(V_ref,6));

disp('-----------pricing results------------');
re = cell(1,nb_dim);
for i=1:11
    re{i} = abs(value{i} - V_ref(i))./V_ref(i);
    disp({'d=' num2str(i+1) vpa(value{i},5) vpa(re{i},3)});
end

%% higher dimensions
% d = 13
value{12} = [2.817529035537446   2.685537188722416];
nb_knots_all{12} = [23297         495       42875]; % ilevel = 4;

% d = 14
value{13} = [2.804341598780551   2.6720];
nb_knots_all{13} = [4117          85       55077;
    30801         561       55077];

% d = 15
value{14} = [2.790995834832607   2.659714602086090];
nb_knots_all{14} = [5021          91       69671;
     40001         631       69671];

% d = 16
value{15} = [2.778358120217755    2.650187960697415];
nb_knots_all{15} = [6049          97       86945;
    51137         705       86945];

re = cell(1,nb_dim);
for i=12:15
    re{i} = abs(value{i} - V_ref(i))./V_ref(i);
    disp({'d=' num2str(i+1) vpa(value{i},5) vpa(re{i},3)});
end