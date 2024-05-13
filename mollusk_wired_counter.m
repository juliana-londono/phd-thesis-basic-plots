%script to create the quadruped gaits + counter sA matrix, where L2 is an
%independent set and new sequence elements are added without adding extra
%nodes to L2, only to counter!

e = 0.07; d = .3;
sequence = [4,6,5,2,3,1,5];
m = size(sequence,2); %number of pulses

%define your layers
%first layer is counter
sAcell{1,1} = [
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0;
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1;
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0;
];
%second layer is independent set
sAcell{2,2} = zeros(6,6);
%third layer is mollusk
sAcell{3,3} = [
    0,0,1,0,0,1;
    1,0,0,1,0,0;
    0,1,0,0,1,0;
    0,0,1,0,0,1;
    1,0,0,1,0,0;
    0,1,0,0,1,0
    ];

%number of layers
N = size(sAcell,2);

%for each layer, get W,b, and X0
for i = 1:N
    Wcell{i,i} = graph2net(sAcell{i,i},e,d);
    n_array(i) = size(Wcell{i,i},1);
end

%layer sizes
n = sum(n_array);

%connections from indep. set to mollusk:
temp = zeros(n_array(3),n_array(2)); %zeros(size of receiver,size of sender)
cols = 1:n_array(2); %sender: indep. set
rows = 1:n_array(3); %receiver: mollusk
for i = 1:size(cols,2)
    temp(rows(i),cols(i)) = 1;
end
sAcell{3,2} = temp;

%connection from counter to indep set:
temp = zeros(n_array(2),n_array(1)); %zeros(size of receiver,size of sender)
cols = [1:2:2*m,2:2:2*m+1]; %sender: counter
rows = [sequence,sequence]; %receiver: indep.set
for i = 1:size(cols,2)
    temp(rows(i),cols(i)) = 1;
end
sAcell{2,1} = temp;

%build the W with feed forward layered structure
for i = 1:N
    for j = 1:N
        if j+1==i %feedforward connections are as specified above
            A = sAcell{i,j};
            Wcell{i,j} = -ones(size(A))+A*e + (A-ones(size(A)))*d;
        elseif i>j | j>i  %everyone else is 0
            Wcell{i,j} = zeros(size(Wcell{i,i},1),size(Wcell{j,j},1));
        end
    end
end

%build the big W from the blocks
for i = 1:size(Wcell,1)
    for j = 1:size(Wcell,1)
        %record the partition indices taus
        if i==j
            tau{i} = (sum(n_array(1:i-1))+1):sum(n_array(1:i));
        end
        W((sum(n_array(1:i-1))+1):sum(n_array(1:i)),(sum(n_array(1:j-1))+1):sum(n_array(1:j))) = Wcell{i,j};
    end
end

figure(1)
imagesc(W)
hold on
for j=1:size(n_array,2)-1
    xline(sum(n_array(1:j))+0.5,'m')
    yline(sum(n_array(1:j))+0.5,'m')
end
hold off
colormap(gray)
colorbar('Ticks',unique(W),...
    'TickLabels',{'-1-\delta','-1+\epsilon','0'})
xticks(1:n)
yticks(1:n)
title(['epsilon = ',num2str(e),', delta = ', num2str(d),...
    ', sequence = ', num2str(sequence)])

n_pulses = m; %as many as length of sequence, to wrap around

counter_nodes = 1:n_array(1);
indep_nodes = (n_array(1)+1):(n_array(1)+n_array(2));
mollusk_nodes = (n_array(1)+n_array(2)+1):(n);

% create pulse vector
%%%%%%%%%%%customizables%%%%%%%%%%%
pulsed_neurons = cell(1,n_pulses);
for i = 1:n_pulses
    pulsed_neurons{i} = 1:(n_array(1)+n_array(2)); %uniform pulse to L1 and L2
end
pulse_separation = 100; %empty so that it is randomly separated by default
w = 3;
h = 5.5;
theta = ones(n,1);
%this initial condition is for exact reproducibility, but mechanism works
%with X0 = 0.1*rand(1,n)
X0 = [0.5,0.5,0.002,0.093,0.069,0.1,0.017,0.014,0.093,0.07,...
    0.007,0.076,0.075,0.092,0.071,0.012,0.002,0.003,0.003,0.025,...
    0.7,0.7,0.7,0.084,0.012,0.028];
% X0 = 0.1*rand(1,n);
X0([counter_nodes(1:2)]) = 0.5;
X0([mollusk_nodes(1:3)]) = 0.7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%silence layer 2 external inputs
theta(tau{2}) = 0;

% [b,T,t,pulse] = make_pulse_train(pulsed_neurons,n,pulse_separation,w,h,theta);
no_pulses = size(pulsed_neurons,2);
T = pulse_separation*ones(1,2*no_pulses+1); %uniformly separated pulses

%pulse times
for i = 2:2:2*no_pulses
    T(i) = w;
end

b = repmat(theta,1,length(T));

for j = 1:no_pulses %for each pulse:
    sigma = pulsed_neurons{j};
    for i = 1:length(sigma) %for each neuron in the j-th group
        b(sigma(i),2*j) = h;
    end
end

%this t,pulse is useful to plot the pulses plot(t,pulse,'-k');
t(1) = 0;
for i=1:size(b,2)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = b(:,i); %doubled theta pulse? yes, to plot that little rectangle
    pulse(:,2*i) = b(:,i);
end
t = t(1:size(pulse,2)); % cut off last value


%make the pulse to L1 smaller!
b(tau{1},2:2:size(b,2)) = 5;

%this t,pulse is useful to plot the pulses plot(t,pulse,'-k');
t(1) = 0;
for i=1:size(b,2)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = b(:,i); %doubled theta pulse? yes, to plot that little rectangle
    pulse(:,2*i) = b(:,i);
end
t = t(1:size(pulse,2)); % cut off last value

f2 = figure(2);
f2.Position = [360 467.6667 560 150.3333];
first_nodes(1) = 1;
for j=2:size(n_array,2)
    first_nodes(j) = sum(n_array(1:j-1))+1;
end
imagesc(b(first_nodes,:))
hold on
for j=1:size(n_array,2)-1
    %xline(sum(n_array(1:j))+0.5,'m')
    yline(j+0.5,'m')
end
hold off
%colormap(gray)
colorbar('Ticks',unique(b))
xticks(1:size(b,2))
xticklabels(arrayfun(@num2str, T, 'UniformOutput', 0))
xlabel('pulse durations')
ylabel('layer')
yticks(1:N)
title('\theta')


%solve ODE
soln = threshlin_ode(W,b,T,X0);

molluskan_colors = [166,206,227;... %baby blue, 1
    178,223,138;...%light green, 2
    251,154,153;...%pink, 3
    31,120,180;...%blue, 4
    51,160,44;...%green, 5
    227,26,28; %red, 6
    ];

%color of indep. set nodes matches the color of aux neuron 1 to 1
indep_nodes_colors = molluskan_colors;

%color of counter also matches
for i = 1:size(sequence,2)
    counter_nodes_colors(i,:) = molluskan_colors(sequence(i),:);
end
%repeat for the other counter nodes the colors
counter_nodes_colors = reshape([counter_nodes_colors(:) counter_nodes_colors(:)]',2*size(counter_nodes_colors,1), []);

colors = [counter_nodes_colors;indep_nodes_colors;molluskan_colors]/255;

% plot solution and pulses
figure(3)
subplot(11,1,1:4)
plot_grayscale(soln.X(:,n:-1:1));
yticklabels(string(20:-5:5))
title(['epsilon = ',num2str(e),', delta = ', num2str(d), '. Sequence: ', num2str(sequence)])
hold on
for j=1:size(n_array,2)-1
    %these divide the layers:
    yline(-sum(n_array(1:j))+n+0.5,'m')
end
hold off

subplot(11,1,5)
plot(t,pulse,'k');
% ax = gca; 
% ax.ColorOrder = colors(1:n,:);
title(['theta pulses with width = ',num2str(w),', height = ', num2str(h)])
xlim([0,max(t)])
ylim([0.5,h+0.5]);
yticks(unique(pulse))

subplot(11,1,6:7)
plot_ratecurves(soln.X(:,mollusk_nodes),soln.time,colors(mollusk_nodes,:));
xlim([0,max(t)])
title('L3 nodes (CPG nodes)')
legend(strsplit(num2str(mollusk_nodes)))

subplot(11,1,8:9)
plot_ratecurves(soln.X(:,indep_nodes),soln.time,colors(indep_nodes,:));
xlim([0,max(t)])
title('L2 nodes')
legend(strsplit(num2str(indep_nodes)))

subplot(11,1,10:11)
plot_ratecurves(soln.X(:,counter_nodes),soln.time,colors(counter_nodes,:));
xlim([0,max(t)])
title('L1 nodes')
legend(strsplit(num2str(counter_nodes)))