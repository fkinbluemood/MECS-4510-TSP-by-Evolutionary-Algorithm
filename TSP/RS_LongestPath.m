% MECS 4510 HOMEWORK1
% Author: Zhengdong Liu Jianfei Pan  UNI:zl2957 jp4201
% This script will import the data of 1000 locations and then it will 
% implement the random search method to find the longest path through all
% points.

% INPUT:   run: number of runs                evl: number of evalutaions 
% OUTPUT:  path_x: x coordinate of path,      path_y:y coordinate of path
%          dx: x coordinate for evaluation    dy: longest distance
%          derr: errorbar

function [path_x, path_y,dx,dy,derr]=RS_LongestPath(runs,evl)

    % import the randomly distributed samples and store them in terms of x and 
    % y coordinates
    Sample=importdata('tsp.txt');
    sample_x=Sample(:,1);
    sample_y=Sample(:,2);

    % Run for 10 times to get the average
    for k=1:runs

        % set the initial path distance for checks 
        dist_i=0;

    % loop over 3e6 iterations to find the shortest path
        for j=1:evl
            x1(j)=j;
            % set the initial distance
            dist=0;
            
            % create an array of the random sequence of the path
            %num=randperm(1000,1000); 
            num=randperm(1000,1000);

            % loop over all points
            for i=1:1000
                if i==1000
                    dist=dist+sqrt( (sample_x(num(1000))-sample_x((num(1))))^2+(sample_y(num(1000))-sample_y((num(1))))^2); 
                else
                    % calculate and add up the total distance
                    dist=dist+sqrt( (sample_x(num(i+1))-sample_x((num(i))))^2+(sample_y(num(i+1))-sample_y((num(i))))^2);
                end
            end

            % update the longest distance 
            if dist>dist_i
                path_order=num;
                dist_final(j)=dist;
                dist_i=dist;
            else
                dist_final(j)=dist_i;
            end
        end
        disp_finalY(:,k)=dist_final;
    end

        new_y=mean(disp_finalY,2);
        sd=std(disp_finalY,[],2);
        err=sd/sqrt(runs);
        dx=linspace(1,evl,10);
        dy=interp1(x1,new_y,dx);
        derr=interp1(x1,err,dx);

    % Loop over points to plot the path
     for i=1:1001
         if i==1001
             path_x(1001)=(sample_x(num(1)));
             path_y(1001)=(sample_y(num(1)));
         else
             path_x(i)=(sample_x(num(i)));
             path_y(i)=(sample_y(num(i)));
         end
     end
end

