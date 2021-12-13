%% modeling PAR-3 segregation by single particle simulation
% main_array arranged by rr (particles0 X tmax (time point) X 3(parameters)
% (parameters: 1 = x, 2 = y, 3 = cluster size n)
% rr is the number of rows, or the number of particles being tracked
% ti is time interval, in one frame per ti seconds
% v_cortical in um/min
%t_totalsegregation in min
function Tracks = MainMatrix_EOD_Dec8th2021 (rr, ti, clusterSize, De, v_cortical, p_off, t_totalsegregation, bins)
%%initialization 
% initialize t and tmax
t = 0;
tmax = fix(t_totalsegregation*60/ti+1);
% p: parameters, 1 = x, 2 = y, 3 = z
parameters = 3; 
Tracks = zeros(rr, tmax, parameters);
%% generate initial cluster position at t=0, assumed randomly distributed.
% generate random x and y positions from -30 to 30
x = rand (rr, 1)*60-30;
Tracks (:,1,1) = x(:,1);
y = rand (rr, 1)*30-15;
Tracks (:,1,2) = y(:,1); 
  for i = 1:rr
    Tracks (i,1,3) = clusterSize;
  end
% display embryo with initalized cluster positions
xi = Tracks (:, 1, 1);
yi = Tracks (:, 1, 2);
ni = Tracks (:, 1, 3);
%figure('Name','Inital Cluster Postion');
%hold on
  %for i = 1:rr
       %if xi(i) > -30 && yi(i)<sqrt((1-(xi(i))^2/30^2)*15^2) && yi(i)> -sqrt((1-(xi(i))^2/30^2)*15^2)
      %scatter (xi (i), yi (i), ni (i)*4, 'filled')
       %end
      % *2 to augment difference in appearance of size
  %end
% draw an embryo in coordinates, assuming that embryo is a 30 by 60um
% eclipse: x^2/25^2+y^2/15^2=1;
a=30; % horizontal radius
b=15; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
%eclipse = plot(x,y,'k')
%drawnow
hold off
%% Segregation
%   figure
%   ah = axes;
   for t = 2 : tmax
      % move every cluster by diffusion (assume same D for now)
      Tracks (:, t, 1:2) = MoveByBrownianDiffusion (Tracks (:, t-1, 1:2), ti, De);
      fprintf('Now analyzing time point %d /%d\n',t, tmax)
      % move every cluster by cortical flow
      for r = 1 : rr
           x_last = Tracks (r, t, 1);
           v = (v_cortical * x_last / 90)+(v_cortical*2 / 3);% linerly degrading cortical flow
             Tracks (r, t, 1) = Tracks (r, t, 1)-(v / 60)*ti;% change gradient vs constrant flow here
             Tracks (r, t, 3) = Tracks (r, t-1, 3);
                          
             %chcek for out of boundary clusters and make appear on the other side
             %if Tracks(r, t, 1) < -30
               %Tracks(r, t, 1)= -30+(-30 - Tracks(r, t, 1));          
                 %else if Tracks(r, t, 1) > 30
                    %Tracks(r, t, 1)= 30-(Tracks(r, t, 1)-30);
                 %end
             %end
             
             %if Tracks(r, t, 2) < -15
               %Tracks(r, t, 2)= 15-(-15 - Tracks(r, t, 2));
                %else if Tracks(r, t, 2) > 15
                %Tracks(r, t, 2)= -15+(Tracks(r, t, 2)-15);
                     %end
             %end
             %if Tracks (r, t, 2) < -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2); 
               % = -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2) - Tracks (r, t, 2);
               %Tracks (r, t, 2) = sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2)-d;
                %if Tracks (r, t, 2) > sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2);
                   %d = Tracks (r, t, 2) - sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2);
                   %Tracks (r, t, 2) = -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2)+d;
                %end
             %end   
      end
      
      %determine which cluster will fall off, and select new binding 
      %to replace these clusters according to density cloud.
      P_allcluster = rand(rr, 1); %generate a random p for each cluster
      index_P = [P_allcluster]< p_off; %compare to p_off to determine whehter drop off membrane
      %choose x cordinates
      distr_x = Tracks(:,t, 1); %fetch current distribution of cluster on x axis
      %calculate cluster size in each bins
      binsize = 60 / bins;
      for b = 1:bins
          loweredge = (b-1)*binsize-30;
          upperedge = b*binsize-30;
          bin_index = [distr_x] > loweredge & [distr_x] < upperedge;
          clusterINbins(b)  = sum(bin_index);
      end
      
      %Assigning new locations
      Xdata = 0:binsize:59.9999999;
      Ydata = clusterINbins;
      for b = 2 : bins-1
          %if Ydata(b) == 0 & Ydata(b+1) == 0
              %Ydata(b)=0;
          %else
          %find the next non-zeros bins
          for j = b+1 : bins-1
              if Ydata(j)>0
                  break;
              end
          end
          for k = b-1 : -1 : 2
              if Ydata(k)>0
                  break;
              end
          end
          Ydata(b) = sum([Ydata(k), Ydata(b), Ydata(j)])/3;
          %end
      end
      Ydata(1) = (Ydata(1)+Ydata(2))/2;
      %Ydata(bins) = (Ydata(bins)+Ydata(bins-1))/2;
      %for b = 2 : bins-1
          %Ydata(b) = sum([Ydata(k), Ydata(b), Ydata(j)])/3;
      %end
      x_prob = Ydata ./ sum(Ydata);
      %disp(x_prob);
      x_whichBIN = randsrc(rr,1,[Xdata; x_prob]);
      x_new = [rand(rr,1)*0.6] + [x_whichBIN]-30;
      y_new = rand (rr, 1)*30-15;
      %figure('Name','x_fit');
      %set(gcf,'Visible','off')
      %x_fit = fitdist(distr_x, 'kernel'); %fitting x coordinates distribution %nonparametric kernel-smoothing distribution
      %index3030 = [x_fit(2).XData] > -30 & [x_fit(2).XData] < 30; %trimming the tails off of fitted curve
      %x_fit(2) is the fitting output from histfit
      %Xdata = x_fit(2).XData(index3030); %Xdata should be bin centers
      %Ydata = x_fit(2).YData(index3030); %Ydata should be the counts in each bin
      %Xdata = -25:1:25; %121 evenly spaced positions along the AP axis
      %Ydata = pdf(x_fit,Xdata);
      %x_prob = Ydata ./ sum(Ydata); %re-calculating the fraction of each bin because the probability do not add up to 1 after trimming the tails
      %x_new = randsrc(rr,1,[Xdata; x_prob]); %generating new x coordinates according to the distribution fit
      %clear x_fit; %tried here if the figure is the problem, NOT helpful
      %generate a random y position according to x
      %y_new = rand (rr, 1)*30-15;
      %for j = 1:rr
          %y_new_boundary = sqrt((1-(x_new(j))^2/30^2)*15^2);
          %y_new(j) = rand(1) * y_new_boundary * 2 - y_new_boundary;
      %end
      Tracks(index_P, t, 1) = x_new(index_P);
      Tracks(index_P, t, 2) = y_new(index_P);
      
%       % display cluster positions
%       cla(ah)
%       hold on;
%       xf = Tracks (:, t, 1);
%       yf = Tracks (:, t, 2);
%       ni = Tracks (:, t, 3);
%       for i = 1:rr 
%           if xf(i) > -30 && yf(i)<sqrt((1-(xf(i))^2/30^2)*15^2) && yf(i)> -sqrt((1-(xf(i))^2/30^2)*15^2)
% 
%           scatter (ah, xf (i), yf(i), ni(i)*4, 'filled')
% 
%           end
%       end
% 
% 
%     % draw an embryo
%     a=30; % horizontal radius
%     b=15; % vertical radius
%     x0=0; % x0,y0 ellipse centre coordinates
%     y0=0;
%     q=-pi:0.01:pi;
%     x=x0+a*cos(q);
%     y=y0+b*sin(q);
%     eclipse = plot(ah, x,y, 'k');
%     drawnow
%     hold off
   end
%% display final cluster positions
  figure('Name','Final Cluster Postion');
  hold on;
  xf = Tracks (:, tmax, 1);
  yf = Tracks (:, tmax, 2);
  ni = Tracks (:, tmax, 3);
  for i = 1:rr 
      if xf(i) > -30 && yf(i)<sqrt((1-(xf(i))^2/30^2)*15^2) && yf(i)> -sqrt((1-(xf(i))^2/30^2)*15^2)
      scatter (xf (i), yf(i), ni (i)*4, 'filled')
      end
  end
  
% draw an embryo
a=30; % horizontal radius
b=15; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
q=-pi:0.01:pi;
x=x0+a*cos(q);
y=y0+b*sin(q);
eclipse = plot(x,y, 'k');
drawnow
hold off
end
%% given MainMatrix(:, t, 1:2), move the clusters according to brownian
function out = MoveByBrownianDiffusion (input, dt, D)
ss = size(input);
%dt = 1; % time 
dx = sqrt(4*D*dt); %displacement in time 
out = input + dx*randn(ss(1), 1, 2);
end

