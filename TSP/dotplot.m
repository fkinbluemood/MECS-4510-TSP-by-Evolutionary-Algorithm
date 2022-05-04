runs=1;
evl=1e4;

[path_x, path_y,dx,dy,derr,x1,new_y]=RS_ShortestPathdot(runs,evl);
 plot(x1,new_y);
 hold on
 scatter(x1,new_y);