% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This function will import the data of 1000 locations and then it will 
% implement the hill climber method to find the longest path through all
% points.

% INPUT:   run: number of runs                evl: number of evalutaions 
% OUTPUT:  path_x: x coordinate of path,      path_y:y coordinate of path
%          dx: x coordinate for evaluation    dy: longest distance
%          derr: errorbar

function [path_x,path_y,dx,dy,derr,dist_final]=HillClimber_longest(runs,evl)
    % import the randomly distributed samples and store them in terms of x and 
    % y coordinates
    Sample=importdata('tsp.txt');
    sample_x=Sample(:,1);
    sample_y=Sample(:,2);
    for k =1: runs
        
        % create the initial order of path 
        num=randperm(1000,1000);
        
        % calculate the intial distance
        dist_i=0;
        for i=1:1000
             if i==1000
                 dist_i=dist_i+sqrt( (sample_x(num(1000))-sample_x((num(1))))^2+(sample_y(num(1000))-sample_y((num(1))))^2);
             else
                 dist_i=dist_i+sqrt( (sample_x(num(i+1))-sample_x((num(i))))^2+(sample_y(num(i+1))-sample_y((num(i))))^2);
             end
        end
        
        % loop over n evaluations to improve the result
        for j=1:evl
            % store the data for x coordinate
            x1(j)=j;
            
            % copy a new order array and randomly swap two adjacent cities
            mutate =num;       
            swapidx=randperm(999,1); %Create random indices for swapping
            mutate(swapidx+1)=num(swapidx); % random swapping
            mutate(swapidx)=num(swapidx+1);
            dist=0;
            % loop over all points, calculate the new distance
            for i=1:1000
                if i==1000
                    dist=dist+sqrt( (sample_x(mutate(1000))-sample_x((mutate(1))))^2+(sample_y(mutate(1000))-sample_y((mutate(1))))^2);
                else
                    dist=dist+sqrt( (sample_x(mutate(i+1))-sample_x((mutate(i))))^2+(sample_y(mutate(i+1))-sample_y((mutate(i))))^2);
                end
            end
            % coompare with initial distance and update the longest order               
            if dist>dist_i
                num=mutate;
                dist_final(j,k)=dist;
                dist_i=dist;
            else
                dist_final(j,k)=dist_i(1);
            end
        end
    end
    
        % calculate the errorbars for these runs
        new_y=mean(dist_final,2);
        sd=std(dist_final,[],2);
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

