close all;
clear;
clc;

sample_point = 1000000;
pd = makedist('Uniform','lower',-pi,'upper',pi);
x = random(pd,sample_point,1);

COS_THETA = cos(x);

%1-a
fc = 2e9;
v = 20*1e3 / (60*60);
c = 3e8;
fm = v*fc / c;

Doppler_Shift = fm .* COS_THETA;

figure(1)
subplot(2, 1, 1);
[num, edge] = histcounts(Doppler_Shift, 'Normalization','pdf', 'BinLimits', [-fm, fm]);
xi = linspace(edge(1), edge(end), length(edge)-1);
plot(xi, num);
title('v=20km/hr fc=2GHz PDF')
xlabel('Doppler Shift (Hz)');
ylabel('Probability Density');

subplot(2, 1, 2);
cdfplot(Doppler_Shift);
title('v=20km/hr fc=2GHz CDF')
xlabel('Doppler Shift (Hz)');
ylabel('Cumulative Probability Density');

%1-b
fc = 26e9;
v = 90*1e3 / (60*60);
c = 3e8;
fm = v*fc / c;

Doppler_Shift = fm * COS_THETA;

figure(2)
subplot(2, 1, 1);
[num, edge] = histcounts(Doppler_Shift, 'Normalization','pdf', 'BinLimits', [-fm, fm]);
xi = linspace(edge(1), edge(end), length(edge)-1);
plot(xi, num);
title('v=90km/hr fc=26GHz PDF')
xlabel('Doppler Shift (Hz)');
ylabel('Probability Density');

subplot(2, 1, 2);
cdfplot(Doppler_Shift);
title('v=90km/hr fc=26GHz CDF')
xlabel('Doppler Shift (Hz)');
ylabel('Cumulative Probability Density');

%1-c
pd_v = makedist('Uniform','lower', 20,'upper',90);
V = random(pd_v,sample_point,1) * 1e3 / (60*60);
fc = 2e9;
FM = V*fc / c;
Doppler_Shift = FM .* COS_THETA;
figure(3)
subplot(2, 1, 1);
[num, edge] = histcounts(Doppler_Shift, 'Normalization','pdf', 'BinLimits', [-max(FM), max(FM)]);
xi = linspace(edge(1), edge(end), length(edge)-1);
plot(xi, num);
title('fc=2GHz PDF')
xlabel('Doppler Shift (Hz)');
ylabel('Probability Density');

subplot(2, 1, 2);
cdfplot(Doppler_Shift);
title('fc=2GHz CDF')
xlabel('Doppler Shift (Hz)');
ylabel('Cumulative Probability Density');

%1-d
fc = 2e9;
v = 20*1e3 / (60*60);
c = 3e8;
fm = v*fc / c;
x = linspace(-1, 1, 1000);
y = 1 ./ (pi*sqrt(1 - (x).^2));
Y = 1 - acos(x) ./ pi;

Doppler_Shift = fm .* COS_THETA;
[num, edge] = histcounts(Doppler_Shift, 'Normalization','pdf', 'BinLimits', [-fm, fm]);
xi = linspace(edge(1), edge(end), length(edge)-1);


figure(4)
subplot(2, 1, 1);
hold on
plot(xi, num);
plot(x*fm, y/fm);
title('v=20km/hr fc=2GHz PDF')
legend('Simulation', 'Theoretical', 'Location', 'northeast');
xlabel('Doppler Shift (Hz)');
ylabel('Probability Density');
hold off
subplot(2, 1, 2);
hold on
cdfplot(Doppler_Shift);
plot(x*fm, Y);
title('v=20km/hr fc=2GHz CDF')
legend('Simulation', 'Theoretical', 'Location', 'southeast');
xlabel('Doppler Shift (Hz)');
ylabel('Cumulative Probability Density');
hold off



fc = 26e9;
v = 90*1e3 / (60*60);
c = 3e8;
fm = v*fc / c;
x = linspace(-1, 1, 1000);
y = 1 ./ (pi*sqrt(1 - x.^2));
Y = 1 - (acos(x)) ./ pi;

Doppler_Shift = fm .* COS_THETA;
[num, edge] = histcounts(Doppler_Shift, 'Normalization','pdf', 'BinLimits', [-fm, fm]);
xi = linspace(edge(1), edge(end), length(edge)-1);

figure(5)
subplot(2, 1, 1);
hold on
plot(xi, num);
plot(x*fm, y/fm);
title('v=90km/hr fc=26GHz PDF')
legend('Simulation', 'Theoretical', 'Location', 'northeast');
xlabel('Doppler Shift (Hz)');
ylabel('Probability Density');
hold off
subplot(2, 1, 2);
hold on
cdfplot(Doppler_Shift);
plot(x*fm, Y);
title('v=90km/hr fc=2GHz CDF')
legend('Simulation', 'Theoretical', 'Location', 'southeast');
xlabel('Doppler Shift (Hz)');
ylabel('Cumulative Probability Density');
hold off
