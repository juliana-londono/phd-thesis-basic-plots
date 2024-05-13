% CTLN counter simulations. For my PhD thesis. 
% written by Carina C and Juliana L

% fancy colors from https://colorbrewer2.org/
colors = [166,206,227;...
31,120,180;...
178,223,138;...
51,160,44;...
251,154,153;...
227,26,28;...
253,191,111;...
255,127,0;...
202,178,214;...
106,61,154;...
177,89,40;...
210,180,140]/255;


% make sA matrices________________________________________________
% graph #1: regular counter, forward edges on odd and even nodes
% graph #2: signed counter, forward edges on odd nodes, backwards on even

sA{1} = [
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0;
];

sA{2} = [
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
];

n = size(sA{1},1);

% make pulse trains__________________________________
p(1) = 3; % pulse width
p(2) = 2;
h(1) = 5; % pulse height
h(2) = 2.5;

% pulse times for graphs 1 & 2
r(1) = 1;
r(2) = 3; % "refractory" period after pulse
for i=1:2
    T{i} = [40 p(i) r(i) 60 p(i) r(i) 80 p(i) r(i) 40 p(i) r(i) ... 
        60 p(i) r(i) 30 p(i) r(i) 50 p(i) r(i) 40];
    nT(i) = length(T{i});
end

for i=1:2
    theta{i} = ones(n,length(T{i}));
end

% pulses for graphs 1 & 2
odds = [1:2:11];
evens = [2:2:12];
theta{1}(:,[2:3:nT(1)]) = h(1);
theta{2}(odds,[2,8,11,14,20]) = h(2);
theta{2}(evens,[5,17]) = h(2);

% during "refractory" period for pulse 2, stimulated neurons get depressed
% this avoids getting stuck in tadpole fixed pts, and stops activity from 
% sliding forward more than one clique
theta{2}(odds,1+[2,8,11,14,20]) = 1-.05;
theta{2}(evens,1+[5,17]) = 1-.05;


% time = [0 T(1) T(1) sum(T(1:2)) sum(T(1:2)) ...
for j=1:2
t(1) = 0;
for i=1:length(T{j})
    t(2*i) = sum(T{j}(1:i));
    t(2*i+1) = sum(T{j}(1:i));
    pulse{j}(2*i-1) = max(theta{j}(:,i));
    pulse{j}(2*i) = max(theta{j}(:,i));
end
t = t(1:length(pulse{j})); % cut off last value
time{j} = t;
end


% solve system and plot results__________________________________
e(1) = .25; d(1) = .5; % for unsigned counter
%e(2) = .51; d(2) = 1.76; % for signed counter
e(2) = .25; d(2) = .5; % for signed counter

% initial conditions on first two neurons 1,2
X0 = zeros(1,n);
X0(1) = .3;
X0(2) = .3;

names = {'CTLN counter','CTLN signed counter'};

for i=1:2
    figure(i)
    % get solution to ode
    soln = sA2soln(sA{i},T{i},X0',e(i),d(i),theta{i});
    subplot(2,1,1)
    plot_grayscale(soln.X);
    title(names{i})
    subplot(2,1,2)
    plot_ratecurves(soln.X,soln.time,colors);
    legend('1','2','3','4','5','6','7','8','9','10','11','12')
    xlim([0,sum(T{i})]);
    subplot(7,1,4)
    plot(time{i},pulse{i},'-k');
    ylim([0,max(pulse{i})+0.5]);
    xlim([0,sum(T{i})]);
    yticks(unique(pulse{i}))
end