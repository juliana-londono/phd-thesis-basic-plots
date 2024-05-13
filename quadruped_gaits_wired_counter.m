%script to create the quadruped gaits + counter sA matrix using the simply
%layered structure (L1 is counter, L2 is cycle, L3 is gaits)

% 13 bound
% 15 pace
% 17 trot
% 19 walk
% 23 pronk

%sequence must be a sequence of auxiliary nodes to stimulate!
sequence = [15,13,23,17,19,23,17]; %this is the one from the paper!
aux_nodes = sort(unique(sequence));
indep_sz = size(aux_nodes,2);
m = size(sequence,2); %number of pulses
e = .25; d = .5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
]; %counter
%second layer is independent set
sAcell{2,2} = zeros(indep_sz,indep_sz);
%third layer is gaits
sAcell{3,3} = [0 1 0 1 1 0 1 0 1 0 0 0 1 0 0 1 1 0 0 0 0 1 0 1;
    1 0 1 0 0 1 0 1 0 1 0 0 0 1 0 1 0 1 1 0 0 0 0 1;
    0 1 0 1 1 0 1 0 0 0 1 0 0 1 1 0 1 0 0 0 1 0 0 1;
    1 0 1 0 0 1 0 1 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 1;
    1 0 1 0 0 1 0 1 1 0 0 0 1 0 0 1 1 0 0 0 0 1 0 0;
    0 1 0 1 1 0 1 0 0 1 0 0 0 1 0 1 0 1 1 0 0 0 0 0;
    1 0 1 0 0 1 0 1 0 0 1 0 0 1 1 0 1 0 0 0 1 0 0 0;
    0 1 0 1 1 0 1 0 0 0 0 1 1 0 1 0 0 1 0 1 0 0 0 0;
    1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0;
    0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0;
    0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 1 1 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 1 1 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0;
    1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0];

%number of layers
N = size(sAcell,2);

% convert to W and get size of layers
for i = 1:N
    Wcell{i,i} = graph2net(sAcell{i,i},e,d);
    n_array(i) = size(Wcell{i,i},1);
end

%network size
n = sum(n_array);

%connections from indep. set to gaits:
temp = zeros(n_array(3),n_array(2)); %zeros(size of receiver,size of sender)
cols = 1:n_array(2); %sender: indep. set
rows = aux_nodes; %receiver: aux nodes of gaits
for i = 1:size(cols,2)
    temp(rows(i),cols(i)) = 1;
end
sAcell{3,2} = temp;

%translate gait sequence to L2 indices:
for i=1:size(sequence,2)
    idx = find(aux_nodes==sequence(i));
    translated_sequence(i) = idx;
end

%connection from counter to indep set:
temp = zeros(n_array(2),n_array(1)); %zeros(size of receiver,size of sender)
cols = [1:2:2*m,2:2:2*m+1]; %sender: counter
rows = [translated_sequence,translated_sequence]; %receiver: indep.set
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

n_pulses = m-1; %one less than number of gaits in sequence, because initial 
% condition will count as one pulse

% create pulse vector
pulsed_neurons = cell(1,n_pulses);
for i = 1:n_pulses
    pulsed_neurons{i} = 1:(n_array(1)+n_array(2)); %uniform pulse to L1 and L2
end
pulse_separation = 48; %empty so that it is randomly separated by default
%w = 2; h = 4.6; %paper parameters??
w = 2; h = 12;
theta = ones(n,1);

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
b(tau{1},2:2:size(b,2)) = round(h/2);

%this t,pulse is useful to plot the pulses plot(t,pulse,'-k');
t(1) = 0;
for i=1:size(b,2)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(:,2*i-1) = b(:,i); %doubled theta pulse? yes, to plot that little rectangle
    pulse(:,2*i) = b(:,i);
end
t = t(1:size(pulse,2)); % cut off last value


% intial condition coincides with first gait in sequence
X0 = zeros(1,n); %all nodes off...
init = 0.3; %...except the ones below initialized at this value
X0([1,n_array(1)+1,n_array(1)+n_array(2)+sequence(1)]) = init;

% X0 = zeros(1,n);
% X0([n_array(1)-1,n_array(1),n_array(1)+n_array(2)+1]) = 0.3;

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
xticks(1:13)
xticklabels(arrayfun(@num2str, T, 'UniformOutput', 0))
xlabel('pulse durations')
ylabel('layer')
yticks(1:N)
title(['b'])


%solve ODE
soln = threshlin_ode(W,b,T,X0);

gait_colors = [202,140,43;...%1
    203,33,40;...%2
    18,128,178;...%3
    42,153,71;...%4
    1,1,1;...%5
    1,1,1;...%6
    1,1,1;...%7
    1,1,1;...%8
    1,1,1;...%9
    1,1,1;...%10
    1,1,1;...%11
    1,1,1;...%12
    241,96,106;...%13
    241,96,106;...%14
    127,197,237;...%15
    127,197,237;...%16
    105,193,134;...%17
    105,193,134;...%18
    255,207,75;...%19
    255,207,75;...%20
    255,207,75;...%21
    255,207,75;...%22
    239,87,37;...%23
    239,87,37]/255;%24

counter_nodes = 1:n_array(1);
indep_nodes = (n_array(1)+1):(n_array(1)+n_array(2));
leg_nodes = (n_array(1)+n_array(2))+1:(n_array(1)+n_array(2))+4;
%gait_specific_aux_nodes = aux_nodes + n_array(1)+n_array(2);
gait_specific_aux_nodes = (n_array(1)+n_array(2))+13:(n_array(1)+n_array(2))+24;


%sequence names
for i = 1:size(sequence,2)
    if sequence(i) == 13
        sequence_names{i} = 'bound';
    elseif sequence(i) == 15
        sequence_names{i} = 'pace';
    elseif sequence(i) == 17
        sequence_names{i} = 'trot';
    elseif sequence(i) == 19
        sequence_names{i} = 'walk';
    elseif sequence(i) == 23
        sequence_names{i} = 'pronk';
    end
end
sequence_names = string(strjoin(sequence_names));

%color of cycle nodes matches the color of CPG neuron
cycle_nodes_colors = gait_colors(aux_nodes,:);

%color of counter nodes matches the color of CPG neuron
for i = 1:size(sequence,2)
    counter_nodes_colors(2*i-1,:) = gait_colors(sequence(i),:);
    counter_nodes_colors(2*i,:) = gait_colors(sequence(i),:);
end

colors = [counter_nodes_colors;cycle_nodes_colors;gait_colors];

% plot solution and pulses
figure(3)
subplot(11,1,1:4)
%plot_grayscale(soln.X(:,n_array(1)+n_array(2)+4:-1:1));
plot_grayscale(soln.X(:,n:-1:1));
%plot_grayscale(soln.X);
hold on
for j=1:size(n_array,2)-1
    %these divide the layers:
    yline(-sum(n_array(1:j))+n+0.5,'m')
    %these divide the gaits:
    %yline(-sum(n_array(1:j))+n+0.5,'m')
end
hold off
%yticklabels(string(n:-5:5))
% title(['epsilon = ',num2str(e),', delta = ', num2str(d), '. Sequence: ', sequence_names])
title(['Sequence: ', sequence_names])

subplot(11,1,5)
plot(t,pulse,'-k');
title(['theta pulses with width = ',num2str(w),', height = ', num2str(h)])
xlim([0,max(t)])
ylim([0.5,h+0.5]);
yticks(unique(pulse))

subplot(11,1,6:7)
plot_ratecurves(soln.X(:,leg_nodes),soln.time,colors(leg_nodes,:));
xlim([0,max(t)])
title(['leg nodes (',num2str(leg_nodes),') of L3'])
legend(strsplit(num2str(leg_nodes)))

% subplot(11,1,8:9)
% plot_ratecurves(soln.X(:,gait_specific_aux_nodes),soln.time,colors(gait_specific_aux_nodes,:));
% xlim([0,max(t)])
% title(['gait-specific auxiliary nodes:', num2str(aux_nodes)])
% legend(strsplit(num2str(gait_specific_aux_nodes)))

subplot(11,1,8:9)
plot_ratecurves(soln.X(:,indep_nodes),soln.time,colors(indep_nodes,:));
xlim([0,max(t)])
title('L2 nodes')
legend(strsplit(num2str(indep_nodes)))

subplot(11,1,10:11)
plot_ratecurves(soln.X(:,[counter_nodes]),...
    soln.time,colors([counter_nodes],:));
xlim([0,max(t)])
title('L1 nodes')
legend(strsplit(num2str([counter_nodes])))