clear all
clc

%% number of iterations (affine DVI)
NIP1_iterNum_mean = [20.050, 21.350, 22.300, 22.550, 23.650, 22.950];
NIP2_iterNum_mean = [27.350, 31.500, 31.000, 47.500, 32.150, 29.944];

IP1_iterNum_mean = [87.150, 140.450, 134.250, 189.000, 187.000, 201.400];
IP2_iterNum_mean = [51.650, 66.050, 81.600, 104.750, 121.900, 142.000];

figure(1)
x = [1, 2, 3, 4, 5, 6];

% iterNum = [NIP1_iterNum_mean', NIP2_iterNum_mean', IP1_iterNum_mean', IP2_iterNum_mean'];
% bar(x, iterNum)

plot(x, NIP1_iterNum_mean, 'LineWidth', 1, 'Marker','*')
hold on
plot(x, NIP2_iterNum_mean, 'LineWidth', 1, 'Marker','^')
hold on
plot(x, IP1_iterNum_mean, 'LineWidth', 1, 'Marker','s')
hold on
plot(x, IP2_iterNum_mean, 'LineWidth', 1, 'Marker','o')

legend({'NIP1', 'NIP2', 'IP1', 'IP2'},'Location','northwest')
xticklabels({'10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'})
xlabel('s^{*}')
ylabel('Number of Iteration');
axis([0.5 6.5 0 210])

%% computation time (affine DVI)
NIP1_time_mean = [0.312, 0.340, 0.346, 0.334, 0.358, 0.346];
NIP2_time_mean = [0.549, 0.664, 0.604, 0.772, 0.574, 0.511];

IP1_time_mean = [0.250, 0.409, 0.392, 0.568, 0.550, 0.596];
IP2_time_mean = [0.135, 0.183, 0.236, 0.304, 0.351, 0.452];


figure(2)
x = [1, 2, 3, 4, 5, 6];

% iterNum = [NIP1_time_mean', NIP2_time_mean', IP1_time_mean', IP2_time_mean'];
% bar(x, iterNum)

plot(x, NIP1_time_mean, 'LineWidth', 1, 'Marker','*')
hold on
plot(x, NIP2_time_mean, 'LineWidth', 1, 'Marker','^')
hold on
plot(x, IP1_time_mean, 'LineWidth', 1, 'Marker','s')
hold on
plot(x, IP2_time_mean, 'LineWidth', 1, 'Marker','o')

legend({'NIP1', 'NIP2', 'IP1', 'IP2'},'Location','northwest')
xticklabels({'10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'})
xlabel('s^{*}')
ylabel('Computation Time [s]');
axis([0.5 6.5 0 0.9])

