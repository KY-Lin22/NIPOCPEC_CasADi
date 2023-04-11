clear all
clc

% number of iterations (cart pole)
NIP1_iterNum_mean = [401.850, 400.400, 452.450, 391.600, 452.600, 407.500];
IP1_iterNum_mean = [1477.500, 1500.550, 1536.350, 1618.350, 1618.900, 1617.750];

% computation time (cart pole)
NIP1_time_mean = [17.101, 18.641, 20.995, 18.298, 20.779, 19.654];
IP1_time_mean = [26.346, 26.892, 26.949, 28.953, 28.618, 29.169];

x = [1, 2, 3, 4, 5, 6];
figure(1)
plot(x, NIP1_iterNum_mean, 'LineWidth', 1, 'Marker','*')
hold on
plot(x, IP1_iterNum_mean, 'LineWidth', 1, 'Marker','s')
legend({'NIP1', 'IP1'},'Location','northwest')
xticklabels({'10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'})
xlabel('s^{*}')
ylabel('Number of Iteration');
axis([0.5 6.5 300 2000])

figure(2)
plot(x, NIP1_time_mean, 'LineWidth', 1, 'Marker','*')
hold on
plot(x, IP1_time_mean, 'LineWidth', 1, 'Marker','s')
legend({'NIP1', 'IP1'},'Location','northwest')
xticklabels({'10^{-3}','10^{-4}','10^{-5}','10^{-6}','10^{-7}','10^{-8}'})
xlabel('s^{*}')
ylabel('Computation Time [s]');
axis([0.5 6.5 15 31])