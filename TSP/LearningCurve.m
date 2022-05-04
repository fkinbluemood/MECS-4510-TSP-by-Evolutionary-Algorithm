% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This script will do :
% 1: plot the shortest path by Ga method (25%_2points)
% 2：Plot the learning curve for shortest distance for following methods:
% random search, hill  climber, beam search, Ga_50%_3points,
% Ga_50%_2points, Ga_25%_2points

% clear workspace and command window
clear;
clc;

% set up the parameters, we run 5 times with a population of 10 and
% iterations of 1e5.
% scheme =2 represents the method finding the longest distance, while
% scheme =1 represents the method finding hte shortest distance
runs=1;
evl=1e5;
population_size=10;
scheme=1;

%plot the learning curve for shortest distance for following methods:
% random search, hill  climber, beam search, Ga_50%_3points,
% Ga_50%_2points, Ga_25%_2points
figure (1)

[~,~,dxRs,dyRs,derrRs]=RS_ShortestPath(runs,evl);
 output_graphHC= errorbar(dxRs,dyRs,derrRs);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#D95319';
hold on

[~,~,dxHc,dyHc,derrHc]=HillClimber(runs,evl);
 output_graphHC= errorbar(dxHc,dyHc,derrHc);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#EDB120';
hold on

[~,~,dxBs,dyBs,derrBs]=BeamSearch(runs,evl);
 output_graphHC= errorbar(dxBs,dyBs,derrBs);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#7E2F8E';

k_point=2;  % set the crossover variation, this time 2 points crossover
 [~, ~,dxEa50,dyEa50,derrEa50]=Ea_50(runs,evl,population_size,scheme,k_point);
output_graphHC= errorbar(dxEa50,dyEa50,derrEa50);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#77AC30';
hold on
 
 k_point=3; % set the crossover variation, this time 2 points crossover
  [~,~,dx2Ea50,dy2Ea50,derr2Ea50]=Ea_50(runs,evl,population_size,scheme,k_point);
 output_graphHC= errorbar(dx2Ea50,dy2Ea50,derr2Ea50);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#0072BD';
hold on
 

[path_xShort,path_yShort,dxEa25,dyEa25,derrEa25]=Ea_25(runs,evl,population_size,scheme);
 output_graphHC= errorbar(dxEa25,dyEa25,derrEa25);
output_graphHC.LineWidth = 1.5;
output_graphHC.Color = '#00FF00';
hold on

xlabel('number of iterations')
ylabel('shortest path distance')
title('Shortest distance Vs No. of iterations')
grid on
 
  legend('Random Search','Hill Climber','Beam Search','Ga50%2Point','Ga50%3Point','Ga25%2Point')

%  figure 2 plots the longest path achieved by Ga_25%_2points
figure (2)
scatter(path_xShort,path_yShort);
hold on
plot(path_xShort,path_yShort,'-o');
grid on
xlabel('x axis')
ylabel('y axis')
title('Shortest path')

