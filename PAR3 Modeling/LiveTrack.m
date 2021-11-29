
% draw an embryo in coordinates, assuming that embryo is a 30 by 50um
% eclipse: x^2/25^2+y^2/15^2=1;
a=30; % horizontal radius
b=15; % vertical radius
x0=0; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
eclipse = plot(x,y)
hold on;


handles = Tracks; %import dataset from MainMatrix
x = handles (:, :, 1); %get x
y = handles (:, :, 2); %get y
n = handles (:, :, 3); %get cluster size
ss = size (x); %get size
N = ss(1);
tt = ss(2);
r=zeros(N,2); % points positions
r(:, 1) = x(:, 1);
r(:, 2) = y(:, 1); %initalize cluster positions
hp=plot(r(1,:),r(2,:),'k.');
%if (get(handles.trackCheckBox,'Value'))
    hold on;
    rhist1 = r(100,:); % position history for particle #1
    rhist2 = r(100,:); % position history for particle #2
    %histp1 = line(rhist1(:,1),rhist1(:,2),'Color','b','LineWidth',1.5);
    %histp2 = line(rhist2(:,1),rhist2(:,2),'Color','r','LineWidth',1.5);
%end
%xlim(handles.axes1,[-L2 L2]);
%ylim(handles.axes1,[-L2 L2]);
%axis(handles.axes1,'equal');
for t = 2:tt
    r (:,1) = x(:, t);
    r (:,2) = y(:, t);
    set(hp,'XData',r(:,1),'YData',r(:,2));
    %if (get(handles.trackCheckBox,'Value'))
        rhist1 = vertcat(rhist1,r(1,:));
        rhist2 = vertcat(rhist2,r(2,:));
        %set(histp1,'XData',rhist1(:,1),'YData',rhist1(:,2));
        %set(histp2,'XData',rhist2(:,1),'YData',rhist2(:,2));
    %end
    xlim([-30,30]);
    ylim([-15,15]);
    pause(0.02)
    
end