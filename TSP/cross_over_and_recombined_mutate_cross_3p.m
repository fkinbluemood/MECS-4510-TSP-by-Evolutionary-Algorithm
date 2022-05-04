% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function implement the crossover operator and recombine priority 
% matrix

% INPUT:   priorities_of_cities:priorities_of_cities  
%          k_point:crossover points
% OUTPUT:  priorities_of_cities: new prioirty matrix


function [priorities_of_cities] = cross_over_and_recombined_mutate_cross_3p(priorities_of_cities,k_point)
%cross_over_and_recombined Summary of this function goes here
% randomly cross over points and then recombine. Crossover can happen among
% two or three parents

number_of_elements = size(priorities_of_cities);
number_of_plans = number_of_elements(1);
priorities_of_cities_crossovered = zeros(number_of_plans,1000);
cross_over_segments = zeros(10,250); % 250 elements for 20 segments

count = 1;
for i = 1:2:2*number_of_plans
    cross_over_segments(i,:) = priorities_of_cities(i-(count-1),251:500);
    cross_over_segments(i+1,:) = priorities_of_cities(i-(count-1),501:750);
    count = count + 1;
end

% randomize the segments and ready to crossover and recombine
cross_over_segments = cross_over_segments(randperm(size(cross_over_segments, 1)), :);

random_cross = randperm(5);

for i = 1:(number_of_plans-2)
    priorities_of_cities(random_cross(i),251:500) = cross_over_segments(2*random_cross(i) -1,:);
    priorities_of_cities(random_cross(i),501:750) = cross_over_segments(2*random_cross(i),:);
end
    %priorities_of_cities(i,:) = priorities_of_cities_crossovered(i,:)/max(priorities_of_cities_crossovered(i,:));
end

