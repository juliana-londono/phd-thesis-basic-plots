% CTLN counter simulations. For my PhD thesis. 
% written by Carina C and Juliana L

% colors = [0 .5 .7; .15 .6 0; .5 .5 .5; .8 .55 0; .8 0 0];
% rearrange: 1 = gold, 2 = red, 3 = blue, 4 = green
% colors = [.8 .55 0; .8 0 0; 0 .5 .7; .15 .6 0];
colors = [];

% sA matrix___________________________________
% patchwork of 6 n=5 limit cycle modules
sA = [
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1;
    1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0;
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0;
];

n = size(sA,1);
e = .51;
d = 1.76;

% pulse times_________________________________________
p = 1; % pulse width
T = [50 p 50 p 50 p 50 p 50 p 50 p 50];
theta = ones(n,length(T));

% initial condition for first 5-star
X0 = zeros(n,1);
X0(1:5) = [.1003; .1002; .1005; .1004; .1009];

%pulses set up
nodes = [6 9 12 15 18 3];
p_times = [2 4 6 8 10 12];

% % uniform stimulation_____________________________________
% % pulses to hit all in-degree 2 nodes: 3, 6, 9, 12, 15, 18
theta(nodes,p_times) = 2.5;

% t = [0 T(1) T(1) sum(T(1:2)) sum(T(1:2)) ...
t(1) = 0;
for i=1:length(T)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = theta(:,i);
    pulse(:,2*i) = theta(:,i);
end
t = t(1:length(pulse)); % cut off last value

% compute solution and plot results______________________________
soln = sA2soln(sA,T,X0,e,d,theta);

figure(1)
subplot(2,1,1)
plot_grayscale(soln.X);
yticks([3,6,9,12,15,18])
title('identical pulses (uniform stimulation)');
subplot(2,1,2)
plot_ratecurves(soln.X,soln.time,colors);
% hold on
% plot(soln.time,sum(soln.X,2),'-k'); % total pop activity
% hold off
xlim([0,sum(T)])
subplot(7,1,4)
plot(t,pulse);
ylim([0,max(pulse,[],'all')+1]);
xlim([0,sum(T)]);
title('pulses hit all in-deg 2 nodes: 3, 6, 9, 12, 15, 18') 

% specific stimulation_____________________________________
theta = ones(n,length(T));
% pulses hit a single node at the time, in the order specifed in array
for i = 1:length(p_times)
    theta(nodes(i),p_times(i)) = 2.5;
end

% t = [0 T(1) T(1) sum(T(1:2)) sum(T(1:2)) ...
t(1) = 0;
for i=1:length(T)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = theta(:,i);
    pulse(:,2*i) = theta(:,i);
end
t = t(1:length(pulse)); % cut off last value

% compute solution and plot results______________________________
soln = sA2soln(sA,T,X0,e,d,theta);

figure(2)
subplot(2,1,1)
plot_grayscale(soln.X);
yticks([3,6,9,12,15,18])
title('targeted pulses (specific stimulation)');
subplot(2,1,2)
plot_ratecurves(soln.X,soln.time,colors);
% hold on
% plot(soln.time,sum(soln.X,2),'-k'); % total pop activity
% hold off
xlim([0,sum(T)])
subplot(7,1,4)
plot(t,pulse);
ylim([0,max(pulse,[],'all')+1]);
xlim([0,sum(T)]);
title('pulses hit a single node at the time') 