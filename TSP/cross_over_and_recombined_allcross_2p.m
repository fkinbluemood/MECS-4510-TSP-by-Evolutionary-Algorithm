function [priorities_of_cities] = cross_over_and_recombined_allcross_2p(priorities_of_cities,k_point)
%cross_over_and_recombined Summary of this function goes here
% randomly cross over points and then recombine. Crossover can happen among
% two or three parents

number_of_elements = size(priorities_of_cities);
number_of_plans = number_of_elements(1);
cross_over_segments = zeros(5,333); % 250 elements for 20 segments

% chopped the inputs into pieces and prepare for the crossover and
% recombination
for i = 1:number_of_plans
    cross_over_segments(i,:) = priorities_of_cities(i,333:665);
end

% randomize the segments and ready to crossover and recombine
cross_over_segments = cross_over_segments(randperm(size(cross_over_segments, 1)), :);

% The next step is to perform cross_over and recombination
for i = 1:number_of_plans
    priorities_of_cities(i,333:665) = cross_over_segments(i,:);
end
end


