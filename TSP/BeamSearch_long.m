% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function will import the data of 1000 locations and then it will 
% implement the parallel-climber method to find the longest path through all
% points.

% INPUT:   run: number of runs                evl: number of evalutaions 
% OUTPUT:  path_x: x coordinate of path,      path_y:y coordinate of path
%          dx: x coordinate for evaluation    dy: longest distance


function [path_x,path_y,dx,dy,derr]=BeamSearch_long(runs,evl)
    % import the randomly distributed samples and store them in terms of x and 
    % y coordinates
    Sample=importdata('tsp.txt');
    sample_x=Sample(:,1);
    sample_y=Sample(:,2);

    % set the initial path distance for checks 
    num_new = zeros(5,1000);
    num = zeros(10,1000);

    % create five random travelling sequence


    for k =1: runs
        % generate the five random path sequences
        for n=1:5
        num(n,:)=randperm(1000,1000);
        end
            % loop over n evaluations to improve the result
            for j=1:evl
                % store the data for x coordinate
                x1(j)=j;
                % Mutate the five sequences and add them into the total population
                % pool. This occurs iteratively. 
                % make a copy of the first five sequences and mutate them afterwards
                mutate =num(1:5,:);
                num(6:10,:)=mutate;
                
                for n=1:5
                    swapidx=randperm(1000,2); %Create random indices for swapping
                    num(n+5,swapidx(1))=num(n,swapidx(2)); % random swapping
                    num(n+5,swapidx(2))=num(n,swapidx(1));
                end

                for m=1:10
                    dist=0;
                    % loop over all points, calculate and add up the total distance,
                    % and store them in dist_final
                    for i=1:1000
                        if i==1000
                            dist=dist+sqrt( (sample_x(num(m,1000))-sample_x((num(m,1))))^2+(sample_y(num(m,1000))-sample_y((num(m,1))))^2);
                        else
                            dist=dist+sqrt( (sample_x(num(m,i+1))-sample_x((num(m,i))))^2+(sample_y(num(m,i+1))-sample_y((num(m,i))))^2);
                        end
                    end
                    % update the longest distance 
                   
                    dist_final(m)=dist;
                end               
                rank_dist=sort(dist_final,'descend'); % return the largest value each row                 
                for n=1:5
                    for l=1:10
                        if rank_dist(n)==dist_final(l)
                            num_new(n,:)=num(l,:); % select top 5 largest dist
                        end
                    end
               
                end
                num(1:5,:)=num_new; % update the sequence
                
                dist_finalNew(j)=rank_dist(1); % store the largest value
            end          
            dist_finalY(:,k)= dist_finalNew; % store thoe largest values  for each run
    end
        % calculate the errorbars for these runs
        new_y=mean(dist_finalY,2);
        sd=std(dist_finalY,[],2);
        err=sd/sqrt(k);
        dx=linspace(1,evl,10);
        dy=interp1(x1,new_y,dx);
        derr=interp1(x1,err,dx);

    % Loop over points to plot the path
     for i=1:1001
         if i==1001
             path_x(1001)=(sample_x(num(1,1)));
             path_y(1001)=(sample_y(num(1,1)));
         else
             path_x(i)=(sample_x(num(1,i)));
             path_y(i)=(sample_y(num(1,i)));
         end
     end

end

