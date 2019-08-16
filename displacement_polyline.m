clc
%clear

axis equal
LO = [];

%Case1
%LD = [0 0; 10 2; 20 0]  %Displaced polyline
%L = [0 3; 9 3; 19 1]    %Polyline

%Case 2
%LD = [3 11; 6 11; 9 11; 9 9; 16 9; 17 12; 20 12];  %Displaced polyline
%L = [10 10; 15 10; 15 13; 10 13; 10 10];     %Polyline
%LD = flip(LD);
% LD =[0 0; 10 0]
% L = [10 10; 15 15]

%Case 3
%LD = [-5 10; 3 11; 6 11; 9 11; 9 9; 16 9; 17 12; 20 12; 25 13; 15 15; 10 15];  %Displaced polyline
%L = [10 10; 15 10; 15 13; 10 13; 10 10];     %Polyline

%Case 4
%LD = [-5 3; -3 3; -2 2; -2 1; -1 1;0 -1; -2 -2; -3 -3; -3 -4; -2 -4; -1 -3; 0 -3; 0 -4; -1 -5; -1 -6; 0 -7;1 -6; 0 -5; 1 -5; 2 -6; 2 -7; 1 -8];
%L = [-3 4; 0 0; 1 -4; 3 -7];

%Case 5
%LD = [-2 8.5; 2 9; 5 9.5; 8 9; 13 9; 17.1 9.3; 20 9; 24 10; 28 9.5];
%L = [10 10; 15 10; 15 13; 10 13; 10 10];

%Case 6
% LD = [-2.2 12; 0 11.9; 2.3 11.5; 4.2 10.7; 6.6 9.3; 9.9 8.4; 12.9 7.4; 16.8 7.1; 19.7 6.9; 20.7 6.9 ];
% L = [6.0 17.3; 4.6 11.4; 14.9 8.5; 17.8 17.2];
% LO = [-2.2 12; 0 11.9; 2.3 11.5; 2.3 10.8; 2.5 10.2; 3.0 9.6; 3.7 9.2; 4.8 8.8; 15.5 5.8; 16.8 7.1; 19.7 6.9; 20.7 6.9] ;

%Case 7
%L2 =[-045.009 -379.789; -003.956 -395.540; -018.916 -441.041; -056.441 -430.509; -045.009 -379.789]
%L = [-72.6653 -431.6708; -65.8289 -433.1175; -56.2316 -446.3816; -70.9090 -457.0015; -84.2346 -440.7453; -72.6653 -431.6708]
%AL = cell(1, 2);
%AL{1} = L;
%AL{2} = L2;
%L =[-045.009 -379.789; -003.956 -395.540; -018.916 -441.041; -056.441 -430.509]
%LD = [2.449 -446.954; -001.530 -447.279; -011.483 -446.487; -026.164 -444.983; -031.905 -443.282; -040.749 -440.387; -049.239 -438.146; -058.118 -433.053; -062.428 -431.460; -066.030 -430.149; -073.282 -428.888; -079.115 -428.642; -081.855 -428.448]

%LD = flip(LD);

%Case 8
LD = [-60.736 -777.082; -76.624 -770.135; -87.882 -766.662; -93.233 -763.835; -101.049 -756.798; -104.116 -750.724; -104.236 -745.311; -102.432 -738.635; -95.939 -731.478; -87.567 -722.531; -77.646 -708.954; -75.000 -702.097; -74.505 -685.768; -73.137 -681.979; -65.561 -671.454;...
     -59.774 -663.455; -57.880 -657.245; -57.775 -649.351; -58.827 -645.562; -61.352 -642.509; -71.453 -634.195; -78.082 -629.984; -81.133 -626.722; -87.972 -617.670; -89.130 -613.144; -90.603 -600.935; -90.287 -590.304; -89.551 -583.778; -88.709 -578.306; -85.552 -570.938; -78.713 -564.728; -68.718 -557.992];
L1 =[-30.840 -776.285; -33.891 -770.707; -36.311 -768.917; -43.255 -765.444; -49.884 -763.970; -77.556 -767.654; -81.821 -767.598; -86.337 -766.537; -92.462 -762.570; -99.785 -755.971; -103.829 -750.649; -103.888 -747.114; -102.475 -743.265; -94.343 -740.795; -85.307 -736.053; ...
     -80.257 -723.508; -74.176 -706.084; -66.995 -685.257; -62.323 -675.771; -58.750 -664.842; -55.623 -657.418; -54.529 -652.023; -64.934 -614.802; -62.872 -608.066; -51.157 -582.256; -42.052 -559.263]
L2 = [-99.555 -585.486; -99.555 -600.396; -94.624 -600.396; -94.624 -585.486; -99.555 -585.486]
L3 = [-79.046 -569.230; -76.780 -572.363; -69.495 -567.094; -71.762 -563.961; -79.046 -569.230]
L4 = [-43.381 -591.612; -38.113 -586.231; -36.011 -588.288; -41.279 -593.670; -43.381 -591.612]
L5 = [-37.332 -580.495; -32.737 -584.531; -29.735 -581.113; -34.330 -577.077; -37.332 -580.495]
L6 = [-45.012 -621.110; -39.633 -625.146; -41.503 -627.639; -46.882 -623.603; -45.012 -621.110]
L7 = [-87.210 -753.316; -81.767 -753.316; -81.767 -749.600; -87.210 -749.600; -87.210 -753.316]
LD2 =[-33.292 -559.311; -37.531 -561.779; -38.071 -566.482; -32.444 -572.418; -28.127 -574.886; -24.119 -579.589; -23.348 -584.137; -27.819 -589.689; -33.215 -594.238; -43.852 -600.945; -51.945 -604.954; -54.258 -613.281; -54.643 -621.685; -53.025 -623.767; -42.927 -631.014; ...
     -38.533 -636.257; -38.456 -639.803; -41.925 -647.436; -48.477 -656.071; -48.400 -661.623; -53.718 -681.360; -64.741 -695.315; -68.749 -704.027; -66.591 -719.447; -73.451 -730.550; -79.541 -736.795; -84.628 -741.035; -97.038 -745.815; -97.655 -752.754; -87.634 -760.183; ...
     -75.521 -762.809; -59.426 -760.606; -54.097 -758.448; -47.842 -753.998; -43.173 -753.954; -34.408 -756.773; -33.747 -763.734; -25.687 -773.702; -23.220 -781.258]    
AL = cell(1, 7);
AL{1} = L1; 
AL{2} = L2; 
AL{3} = L3; 
AL{4} = L4; 
AL{5} = L5; 
AL{6} = L6; 
AL{7} = L7; 
LD = flip(LD);

%Input parameters
dmax = 20;
dbuff = 6.475;
ten = 0.0;
smooth = 0;
ds = 0;

%Densify polyline
[LDs, IDs] = densifyPolyline(LD, ds);
LDO = LD;
LD = LDs;

hold on

%Select method
methods = {'M1: normal, '; 'M2: average bisector, '; 'M3: average weighted bisector, '; 'M4: weighted bisector, greedy, '};
method = 4; func = 1;

%New polyline: displacement
if (method == 1)
    %LN = displaceByMinDist(L, LD, dmax, dbuff, func);
    %LN = displaceByMinDist2(L, LD, IDs, dmax, dbuff, ten, smooth);
    LN = displaceByMinDist2AL(AL, LD, IDs, dmax, dbuff, ten, smooth);
%elseif (method == 2)
%    LN = displaceByMinDistAndBis(L, LD, dmax, dbuff, func);
%elseif (method == 3)
%    LN = displaceByMinDistAndWeightedBis(L, LD, dmax, dbuff, func);
elseif (method == 4)
    %LN = displaceByWeightedBisGreedy(L, LD, IDs, dmax, dbuff, ten, func);
    %LN = displaceByWeightedBisGreedy2(L, LD, IDs, dmax, dbuff, ten, smooth);
    %LN = displaceByWeightedBisGreedy3(L, LD, IDs, dmax, dbuff, ten, smooth);
    %LN = displaceByWeightedBisGreedy3AL(AL, LD, IDs, dmax, dbuff, ten, smooth);
    LN = displaceByWeightedBisGreedy4AL(AL, LD, IDs, dmax, dbuff, ten, smooth);
end

% if (length(AL) > 0)
    %Plot buffer
    buff1 = []; buff2 = [];
    for i = 1 : length(AL)
        buff1 = polybuffer(AL{i},'lines', dbuff);
        buff2 = polybuffer(AL{i},'lines', dmax);

        plot(buff2, 'FaceColor', 'y', 'EdgeColor', 'y', 'FaceAlpha',0.3, 'DisplayName','Buffer outer, propagation');
        plot(buff1, 'FaceColor', [0.9290 0.6940 0.1250], 'EdgeColor', [0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'DisplayName','Buffer inner, displacement');
    end

    %Plot original polyline and barriers
    h2 = plot (LDO(:, 1), LDO(:, 2), 'b-*', 'LineWidth', 2, 'DisplayName','Source polyline');

    for i = 1 : length(AL)
        plot(AL{i}(:, 1), AL{i}(:, 2), 'k', 'LineWidth', 2, 'DisplayName','Barrier');
    end

    %Plot new line
    plot(LN(:, 1), LN(:, 2), 'r-*', 'LineWidth', 2, 'DisplayName','Displaced polyline');

    %Plot old solution
    if (length(LO) > 1) 
        plot(LO(:, 1), LO(:, 2), 'g', 'LineWidth', 2, 'DisplayName','Old solution');
    end;
    
%  else
%     %Plot buffer
%     buff1 = polybuffer(L,'lines', dbuff);
%     buff2 = polybuffer(L,'lines', dmax);
%     h1 = plot(buff2, 'FaceColor', 'y', 'EdgeColor', 'y', 'FaceAlpha',0.3, 'DisplayName','Buffer outer, propagation');
%     h0 = plot(buff1, 'FaceColor', [0.9290 0.6940 0.1250], 'EdgeColor', [0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'DisplayName','Buffer inner, displacement');
% 
%     %Plot original polyline and barrier
%     h2 = plot (LDO(:, 1), LDO(:, 2), 'b-*', 'LineWidth', 2, 'DisplayName','Source polyline');
%     h3 = plot(L(:, 1), L(:, 2), 'k', 'LineWidth', 2, 'DisplayName','Barrier');
% 
%     %Plot new line
%     h4 = plot(LN(:, 1), LN(:, 2), 'r-*', 'LineWidth', 2, 'DisplayName','Displaced polyline');
% 
%     %Plot old solution
%     if (length(LO) > 1) 
%         h5 = plot(LO(:, 1), LO(:, 2), 'g', 'LineWidth', 2, 'DisplayName','Old solution');
%     end;
% 
%     axis equal
%     title_text = strcat(methods(method, :), strcat('buff. in= ', strcat( num2str(dbuff), strcat(' m, ', 'buff. out= ', strcat(num2str(dmax), 'm', strcat(', tension= ', num2str(ten)))))));
%     title(title_text);
%     legend([h0, h1, h2, h3, h4]);
%     %legend([h0, h1, h2, h3]);    
% end
   
axis equal
%title_text = strcat(methods(method, :), strcat('buff. in= ', strcat( num2str(dbuff), strcat(' m, ', 'buff. out= ', strcat(num2str(dmax), 'm', strcat(', tension= ', num2str(ten)))))));
%title(title_text);
%legend([h0, h1, h2, h3, h4]);
%legend([h0, h1, h2, h3]);