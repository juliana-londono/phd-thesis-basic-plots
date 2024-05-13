% octahedral network for molluskan hunting directions (Clione)

sA = [0,0,1,0,0,1;1,0,0,1,0,0;0,1,0,0,1,0;0,0,1,0,0,1;1,0,0,1,0,0;0,1,0,0,1,0];
n = size(sA,1);

colors = (1/255)*[166,206,227;... %baby blue, 1
    178,223,138;...%light green, 2
    251,154,153;...%pink, 3
    31,120,180;...%blue, 4
    51,160,44;...%green, 5
    227,26,28]; %red, 6

%values to play around with%
e = .25; d = .5;
%e = 0.4; d = 0.75;
e = 0.07; d = .3; %wired counter parameters

X0 = 0.01*rand(n,1);
X0([1,5,6]) = 0.2;

no_pulses = 5;
pulsed_neurons = [4,2,5,3,6];
sigmas = num2cell(pulsed_neurons'); %this is unnecesary, fix it!
p = 2; % pulse width
a = 1.25; %pulse height
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = round(100*rand(1,2*no_pulses+1))+50; %randomly separated pulses
for i = 2:2:2*no_pulses
    T(i) = p;
end

theta = 1*ones(n,length(T));

for j = 1:no_pulses %for each pulse:
    sigma = sigmas{j};
    for i = 1:length(sigma) %for each neuron in the j-th group
        theta(sigma(i),2*j) = a;
    end
end

% pulse times
t(1) = 0;
for i=1:size(theta,2)
    t(2*i) = sum(T(1:i));
    t(2*i+1) = sum(T(1:i));
    pulse(2*i-1) = max(theta(:,i)); %doubled theta pulse? yes to plot that little rectangle
    pulse(2*i) = max(theta(:,i));
end
t = t(1:length(pulse)); % cut off last value

for i = 1:length(pulsed_neurons)
    theta_pulse_colors(i,:) = colors(pulsed_neurons(i),:);
    greyscale_pulse_colors(2*i-1,:) = colors(pulsed_neurons(i),:);
    greyscale_pulse_colors(2*i,:) = colors(pulsed_neurons(i),:);
end
pulse_lines_greyscale = unique(100*t(2:length(t)-1)); %has beginning and end of pulse
pulse_lines_ratecurves = unique(t(2:length(t)-1));

% get solution to ode
soln = sA2soln(sA,T,X0,e,d,theta);

figure(1)
%%%%%
subplot(5,1,1:2)
plot_grayscale(soln.X);
hold on
for i = 1:2*no_pulses %two lines per pulse: beggining and end
xline(pulse_lines_greyscale(i),'Color',greyscale_pulse_colors(i,:),'LineWidth',1,'Alpha',1) %plot red line on greyscale where there is a pulse
end
hold off
yticks(1:n)
yticklabels({'1.up','2.left','3.front','4.down','5.right','6.back'})
title(['molluskan hunting directions with \epsilon = ',num2str(e),', \delta = ', num2str(d)])
%%%%%
subplot(5,1,3)
plot(t,pulse,'-k');
title(['theta pulses with width = ',num2str(p),', height = ', num2str(a)])
xlim([0,max(t)])
ylim([min(pulse)-0.1,max(pulse)+0.1]);
v = [t',pulse'];
v = v(2:length(t)-1,:);
f = 1:4*no_pulses;
f = reshape(f,4,no_pulses)';
patch('Faces',f,'Vertices',v,'FaceVertexCData',theta_pulse_colors,'FaceColor','flat')
%%%%
subplot(5,1,4:5)
plot_ratecurves(soln.X,soln.time,colors);
hold on
for i = 1:no_pulses %plot only beggining of pulse
xline(pulse_lines_ratecurves(2*i-1),'Color',greyscale_pulse_colors(2*i-1,:),'LineWidth',1,'Alpha',1,'LineStyle',':') %plot red line on ratecurves where there is a pulse
xline(pulse_lines_ratecurves(2*i),'Color',greyscale_pulse_colors(2*i,:),'LineWidth',1,'Alpha',1,'LineStyle',':') %plot red line on ratecurves where there is a pulse
antagonist = rem(pulsed_neurons(i)+3,6); %opposite neuron is three neurons more mod 6
if antagonist == 0 % 6 mod 6 is 0, so must change it to 6
    antagonist = 6;
end
%save end and beggining of pulse values:
[~,soln_idx] = min(abs(soln.time-pulse_lines_ratecurves(2*i-1))); %convert to soln.X index 
beginning_of_pulse_values(i,:) = soln.X(soln_idx,:);
%
[~,soln_idx] = min(abs(soln.time-pulse_lines_ratecurves(2*i))); %convert to soln.X index 
end_of_pulse_values(i,:) = soln.X(soln_idx,:);
end
hold off
xlim([0,max(t)])
legend('1.up','2.left','3.front','4.down','5.right','6.back')
%ylim([0,0.5]);

