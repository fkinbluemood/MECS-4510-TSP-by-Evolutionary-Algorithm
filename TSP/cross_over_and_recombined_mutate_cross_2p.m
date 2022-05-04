% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function implement the crossover operator and recombine priority 
% matrix

% INPUT:   priorities_of_cities:priorities_of_cities  
%          k_point:crossover points
% OUTPUT:  priorities_of_cities: new prioirty matrix

function [priorities_of_cities] = cross_over_and_recombined_mutate_cross_2p(priorities_of_cities, k_point)
%cross_over_and_recombined Summary of this function goes here
% randomly cross over points and then recombine. Crossover can happen among
% two or three parents

number_of_elements = size(priorities_of_cities);
number_of_plans = number_of_elements(1);
cross_over_segments = zeros(number_of_plans,333); % 250 elements for 20 segments

% chopped the inputs into pieces and prepare for the crossover and
% recombination
for i = 1:number_of_plans
    cross_over_segments(i,:) = priorities_of_cities(i,333:665);
end

% randomize the segments and ready to crossover and recombine
cross_over_segments = cross_over_segments(randperm(size(cross_over_segments, 1)), :);

% The next step is to perform cross_over and recombination
random_cross = randperm(number_of_plans);

for i = 1:(number_of_plans*0.7)
    priorities_of_cities(random_cross(i),333:665) = cross_over_segments(random_cross(i),:);
end
end


