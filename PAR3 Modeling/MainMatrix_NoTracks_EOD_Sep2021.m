

%% modeling PAR-3 segregation by single particle simulation
% main_array arranged by rr (particles0 X tmax (time point) X 3(parameters)
% (parameters: 1 = x, 2 = y, 3 = cluster size n)

% rr is the number of rows, or the number of particles being tracked
% ti is time interval, in one frame per ti seconds
% v_cortical in um/min

function Tracks = MainMatrix_NoTracks_EOD_Sep2021 (rr, ti, clusterSize, De, v_cortical);
%initializing t and tmax
t = 0;
tmax = fix(7*60/ti+1);
% p: parameters, 1 = x, 2 = y, 3 = z
parameters = 3; 
Tracks = zeros(rr, tmax, parameters);


 
%% initialize cluster position at t=0, assumed randomly distributed.
% generate random x and y positions from -30 to 30
x = rand (rr, 1)*60-30;
Tracks (:,1,1) = x(:,1);
% loop to generage random y position for each track
i = 1;
  for i = 1:rr
  y_range_positive = sqrt((1-(x(i, 1))^2/30^2)*15^2);
  Tracks ( i, 1, 2) = rand (1, 1) * y_range_positive * 2 - y_range_positive;
  end
  
%% generate a random inital cluster size, assuming random and unchanged with t

for i = 1:rr
    %n = randsample (10, 1);
    Tracks (i,1,3) = clusterSize;
end

%% display embryo with initalized cluster positions
xi = Tracks (:, 1, 1);
yi = Tracks (:, 1, 2);
ni = Tracks (:, 1, 3);
figure('Name','Inital Cluster Postion');
hold on
  for i = 1:rr
      scatter (xi (i), yi (i), ni (i)*4, 'filled')
      % *2 to augment difference in appearance of size
  end


% draw an embryo in coordinates, assuming that embryo is a 30 by 60um
% eclipse: x^2/25^2+y^2/15^2=1;
a=30; % horizontal radius
b=15; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
eclipse = plot(x,y)
hold off
  %% Segregation
  for t = 2 : tmax
      % move every cluster by diffusion (assume same D for now)
      Tracks (:, t, 1:2) = MoveByBrownianDiffusion (Tracks (:, t-1, 1:2), ti, De);
      % new version, no threshold, each cluster size will have its own De
      for r = 1 : rr
           x_last = Tracks (r, t, 1);
           v = (v_cortical * x_last / 60)+(v_cortical / 2);
             Tracks (r, t, 1) = Tracks (r, t, 1)-(v / 60)*ti;% change gradient vs constrant flow here
             Tracks (r, t, 3) = Tracks (r, t-1, 3);
          
             %chcek out of boundary and make appear on the other side
             if Tracks(r, t, 1) < -30
               Tracks(r, t, 1)= -30+(-30 - Tracks(r, t, 1));
             
           else if Tracks(r, t, 1) > 30
                    Tracks(r, t, 1)= 30-(Tracks(r, t, 1)-30);
               end
           end
           if Tracks(r, t, 2) < -15
               Tracks(r, t, 2)= 15-(-15 - Tracks(r, t, 2));
           else if Tracks(r, t, 2) > 15
               Tracks(r, t, 2)= -15+(Tracks(r, t, 2)-15);
               end
           end

          if Tracks (r, t, 2) < -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2); 
               d = -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2) - Tracks (r, t, 2);
               Tracks (r, t, 2) = sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2)-d;
              if Tracks (r, t, 2) > sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2);
                   d = Tracks (r, t, 2) - sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2);
                   Tracks (r, t, 2) = -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2)+d;
               end
          end
      end
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
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
eclipse = plot(x,y, 'k')
hold off

%clean up tacks that are out of boundaries
%for r = 1:rr
%    for t = 2:tmax
%       if Tracks(r, t, 1) < -30 || Tracks (r, t, 2) < -sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2) || Tracks (r, t, 2) > sqrt((1-(Tracks(r, t,1))^2/30^2)*15^2);
%          Tracks (r, t, :) = NaN;
%        end
%    end
%end 

%calculate percentage of clusters in the anterior domain
%Ant = 0;
%Post = 0;
%for r = 1:rr
%    if Tracks (r, tmax, 1) < 0;
%        Ant =  Ant +1;
%    else if Tracks (r, tmax, 1) > 0;
%            Post = Post +1;
%        percentAnt = Ant./(Ant + Post)    
%   end
%    end
%end

%InitialAnt = 0;
%InitialPost = 0;
%for r = 1:rr
%    if Tracks (r, 1, 1) < 4.2;
%        InitialAnt =  InitialAnt +1;
%    else if Tracks (r, 1, 1) > 4.2;
%            InitialPost = InitialPost +1;
%        percentAntinital = InitialAnt./(InitialAnt + InitialPost)    
%    end
%    end
%end


end

%% given MainMatrix(:, t, 1:2), move the clusters according to brownian
function out = MoveByBrownianDiffusion (input, dt, D)
ss = size(input);
%dt = 1; % time 
dx = sqrt(4*D*dt); %displacement in time 
out = input + dx*randn(ss(1), 1, 2);
end