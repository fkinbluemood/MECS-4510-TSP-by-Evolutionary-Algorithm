% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function get the lastest travel plan based on priority of cities 

% INPUT:   priorities_of_cities:priorities_of_cities  

% OUTPUT:  traversing_sequence: travel plan


function [traversing_sequence] = get_travelPlan(priorities_of_cities)
% get_travelPlan Summary of this function goes here
% Detailed explanation goes here
number_of_elements = size(priorities_of_cities);
number_of_plans = number_of_elements(1);
traversing_sequence = zeros(number_of_plans,1000);

for i = 1:number_of_plans
    
[~,traversing_sequence(i,:)] = sort(priorities_of_cities(i,:),'descend');

end

