% import the randomly distributed samples and store them in terms of x and 
% y coordinates
clear;
clc;
Sample=importdata('tsp.txt');
size(Sample)
sample_x=Sample(:,1);
sample_y=Sample(:,2);
%cluster points and them recombine them, use kmeans clustering with k = 4.
group1_x = [];
group1_y = [];
group2_x = [];
group2_y = [];
group3_x = [];
group3_y = [];
group4_x = [];
group4_y = [];

[idx, C] = kmeans(Sample,4);
for i = 1:1000
    if idx(i) == 1
        group1_x = [group1_x,sample_x(i)];
        group1_y = [group1_y,sample_y(i)];
    elseif idx(i) == 2
        group2_x = [group2_x, sample_x(i)];
        group2_y = [group2_y, sample_y(i)];
    elseif idx(i) == 3
        group3_x = [group3_x, sample_x(i)];
        group3_y = [group3_y, sample_y(i)];
    else
        group4_x = [group4_x, sample_x(i)];
        group4_y = [group4_y, sample_y(i)];
    end
end

sample_x = [group1_x, group2_x, group3_x, group4_x];
sample_y = [group1_y, group2_y, group3_y, group4_y];
size_of_cities = [size(group1_x),size(group2_x),size(group3_x),size(group4_x)]
k1 = size(group1_x,2);
k2 = size(group1_x,2) + size(group2_x,2)+ size(group3_x,2);

% Set some hyperparameters
runs=1;
evl=1e6;
population_size = 20;

% set the initial path distance for checks 
num_new = zeros(population_size/2,1000);
num = zeros(population_size,1000);
priorities_of_cities = zeros(population_size,1000);

for k =1: runs
    
    % NEW: Add weights to all cities and initialize them with random
    % priorities and normalize them. Each priority value corresponds to the
    % city in the smae canonical order. High priority means the particular city
    % will be visited earlier than the city with low priority.
    
    for n=1:population_size/2
    priorities_of_cities(n,:) = [750+randperm(250,size(group1_x,2)),500+randperm(250,size(group2_x,2)), 250+randperm(250,size(group3_x,2)), randperm(250,size(group4_x,2))];
%     priorities_of_cities(n,:) = randperm(1000,1000);
    priorities_of_cities(n,:) = priorities_of_cities(n,:)/max(priorities_of_cities(n,:)) ;% Normalize
    %num(n,:)=randperm(1000,1000);
    end

        % loop over n evaluations to improve the result
        for j=1:evl
            % store the data for x coordinate
            x1(j)=j;
            % Recombination
            mutate =priorities_of_cities(1:population_size/2,:);
            priorities_of_cities(population_size/2 +1:population_size,:)=mutate;
            % Mutation Starts
            % Mutate the five sequences and add them into the total population
            % pool. This occurs iteratively. 
            % make a copy of the first five sequences and mutate them afterwards
            for n=1:population_size/2
                swapidx=randperm(1000,2); %Create random indices for swapping
%                 swapidx_1 = randperm(250,2);
% %               swapidx_4 = 750 + randperm(250,2);
%                 priorities_of_cities(n+population_size/2,swapidx_1(1))=priorities_of_cities(n,swapidx_1(2)); % random swapping
%                 priorities_of_cities(n+population_size/2,swapidx_1(2))=priorities_of_cities(n,swapidx_1(1));
%                 priorities_of_cities(n+population_size/2,swapidx_4(1))=priorities_of_cities(n,swapidx_4(2)); % random swapping
%                 priorities_of_cities(n+population_size/2,swapidx_4(2))=priorities_of_cities(n,swapidx_4(1));
                priorities_of_cities(n+population_size/2,swapidx(1))=priorities_of_cities(n,swapidx(2)); % random swapping
                priorities_of_cities(n+population_size/2,swapidx(2))=priorities_of_cities(n,swapidx(1));
            end
            %Mutation Ends
            %Crossover starts, here we use three-point crossover on the priorities to generate new travel plans. Using 
            %three-point is because we identify that the cities are from four major regions.

            k_point = 2;
            [priorities_of_cities(population_size/2 + 1 :population_size,:)] = cross_over_and_recombined_mutate_cross_2pk(priorities_of_cities(population_size/2 + 1 :population_size,:),k1,k2);

            % Find the travel plan (an 10*1000 array of indices) based on
            % the priorities, get_travelPlan returns the indices that
            % determine the traversing sequence
            [travel_plan] = get_travelPlan(priorities_of_cities);
            
            for m=1:population_size
                dist=0;
                % loop over all points, calculate and add up the total distance,
                % and store them in dist_final
                for i=1:1000
                    if i==1000
                        dist=dist+sqrt( (sample_x(travel_plan(m,1000))-sample_x((travel_plan(m,1))))^2+(sample_y(travel_plan(m,1000))-sample_y((travel_plan(m,1))))^2);
                    else
                        dist=dist+sqrt( (sample_x(travel_plan(m,i+1))-sample_x((travel_plan(m,i))))^2+(sample_y(travel_plan(m,i+1))-sample_y((travel_plan(m,i))))^2);
                    end
                end
                % update the shortest distance 

                dist_final(m)=dist;
            end   
            
            [rank_dist,rank_dist_index]=sort(dist_final); % return the smallest value each row 
            %priorities_of_cities
            for n = 1:population_size/2 + 1
                priorities_of_cities_selected(n,:)=priorities_of_cities(rank_dist_index(n),:); % select top 5 smallest dist
            end
            %priorities_of_cities
            priorities_of_cities(1:(population_size/2 + 1),:)=priorities_of_cities_selected; % update the sequence
            dist_finalNew(j)=rank_dist(1); % store the shortest value
        end          
        dist_finalY(k,:)= dist_finalNew; % store thoe shortesst values  for each run
end
%     % calculate the errorbars for these runs
%     new_y=mean(dist_finalY,2);
%     sd=std(dist_finalY,[],2);
%     err=sd/sqrt(k);
%     dx=linspace(1,evl,10);
%     dy=interp1(x1,new_y,dx);
%     derr=interp1(x1,err,dx);

% Loop over points to plot the path
 for i=1:1001
     if i==1001
         path_x(1001)=(sample_x(travel_plan(1,1)));
         path_y(1001)=(sample_y(travel_plan(1,1)));
     else
         path_x(i)=(sample_x(travel_plan(1,i)));
         path_y(i)=(sample_y(travel_plan(1,i)));
     end
 end

% plot the shortest path achieved by Parallel Climber
figure(1)
plot(path_x,path_y,'-o');
xlabel('x-axis')
ylabel('y-axis')
title('The Shortest Path by the Genetic Algorithm')
