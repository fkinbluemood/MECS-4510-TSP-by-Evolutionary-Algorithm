% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function will plot the bar graph comparing the results obtained by 
% using Random search and Ga_25%_2points

% clear command window and workspace
clear;
clc;

% import the collected data 
  vals = [86.031 463.591 ; 789.498 539.855 ];
x = categorical({'Shortest distance','Longest distance'});
x = reordercats(x,{'Shortest distance','Longest distance'});
b = bar(x,vals);




xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData); 
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
legend('Ga method','Random Search')
grid on
ylabel('Total distance')