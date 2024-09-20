close all;
clear;
clc;

Es_N0 = [1, 3, 5, 7, 9];    %Symbol Energy to Noise Ratio
L = [1, 2, 3, 4];   %Diversity Branch

N = 1e5;  %Number of symbols
Es = 1; %Symbol Energy
es_n0 = 10.^(Es_N0/10);
N0 = 1./es_n0; %Noise Energy

BER_SC = zeros(length(L), length(Es_N0));
BER_MRC = zeros(length(L), length(Es_N0));
BER_EGC = zeros(length(L), length(Es_N0));
BER_DC = zeros(length(L), length(Es_N0));

%Rayleigh fading parameters
K = 0;
mean = sqrt(K/(K+1));
sigma = sqrt(1/(2*(K+1)));

%Generate QPSK symbols
data = rand(2,N);
data = (2*(data > 0.5)-1) / sqrt(2);


for i = 1:length(Es_N0)
    
    %r = g*s + n
    for j = 1 : length(L)
        %Generate Noise & Fading Gain
        n = normrnd(0, sqrt(N0(i)/2), 2, N, L(j)) + 1i*normrnd(0, sqrt(N0(i)/2), 2, N, L(j));
        
        g = normrnd(mean, sigma, 1, N, L(j)) + 1i*normrnd(mean, sigma, 1, N, L(j));
        g = repmat(g,2,1,1);
        
        s = repmat(data, 1, 1, L(j));
        r = g.*s + n;
    
    
        %Selective Combining
        g_magnitude = abs(g(1,:,:));
        [~, selected_branch] = max(g_magnitude, [], 3);
        for k = 1:N
            g_sc(:,k) = g(:,k,selected_branch(k));
            r_sc(:,k) = r(:,k,selected_branch(k));
        end
        r_sc = r_sc .* exp(-1i*angle(g_sc));


        %Maximal Ratio Combining
        r_mrc = sum(conj(g).*r, 3);

        %Equal Gain Combining
        r_egc = sum(exp(-1i*angle(g)).*r,3);

        %Direct Combining
        r_dc = sum(r,3);
        g_dc = sum(g,3);
        r_dc = r_dc .* exp(-1i*angle(g_dc));


        %Demodulation
        data_sc = sign(real(r_sc))./sqrt(2);
        data_mrc = sign(real(r_mrc))./sqrt(2);
        data_egc = sign(real(r_egc))./sqrt(2);
        data_dc = sign(real(r_dc))./sqrt(2);

        %Calculate BER
        BER_SC(j, i) = sum(data_sc(:) ~= data(:)) / N;
        BER_MRC(j, i) = sum(data_mrc(:) ~= data(:)) / N;
        BER_EGC(j, i) = sum(data_egc(:) ~= data(:)) / N;
        BER_DC(j, i) = sum(data_dc(:) ~= data(:)) / N;
    
    end
    
end

%Plot BER versus Es/N0 curve (Rayleigh)
figure
semilogy(Es_N0, BER_SC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_SC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_SC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_SC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Selective Combining (Rayleigh)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_MRC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_MRC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_MRC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_MRC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Maximal Ratio Combining (Rayleigh)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_EGC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_EGC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_EGC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_EGC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Equal Gain Combining (Rayleigh)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_DC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_DC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_DC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_DC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Direct Combining (Rayleigh)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

%Comparison among all combining strategies
figure
semilogy(Es_N0, BER_SC(2,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_MRC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_EGC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_DC(2,:),'-o','linewidth',1);
legend('Selective Combining','Maximal Ratio Combining','Equal Gain Combining','Direct Combining');
title('BER for QPSK with SC & MRC & EGC & DC (Rayleigh)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

BER_MRC_Rayleigh = BER_MRC;

%Ricean fading parameters
K = 1;
mean = sqrt(K/(K+1));
sigma = sqrt(1/(2*(K+1)));

for i = 1:length(Es_N0)
    
    %r = g*s + n
    for j = 1 : length(L)
        %Generate Noise & Fading Gain
        n = normrnd(0, sqrt(N0(i)/2), 2, N, L(j)) + 1i*normrnd(0, sqrt(N0(i)/2), 2, N, L(j));
        
        g = normrnd(mean, sigma, 1, N, L(j)) + 1i*normrnd(mean, sigma, 1, N, L(j));
        g = repmat(g,2,1,1);
        
        s = repmat(data, 1, 1, L(j));
        r = g.*s + n;
    
    
        %Selective Combining
        g_magnitude = abs(g(1,:,:));
        [~, selected_branch] = max(g_magnitude, [], 3);
        for k = 1:N
            g_sc(:,k) = g(:,k,selected_branch(k));
            r_sc(:,k) = r(:,k,selected_branch(k));
        end
        r_sc = r_sc .* exp(-1i*angle(g_sc));


        %Maximal Ratio Combining
        r_mrc = sum(conj(g).*r, 3);

        %Equal Gain Combining
        r_egc = sum(exp(-1i*angle(g)).*r,3);

        %Direct Combining
        r_dc = sum(r,3);
        g_dc = sum(g,3);
        r_dc = r_dc .* exp(-1i*angle(g_dc));


        %Demodulation
        data_sc = sign(real(r_sc))./sqrt(2);
        data_mrc = sign(real(r_mrc))./sqrt(2);
        data_egc = sign(real(r_egc))./sqrt(2);
        data_dc = sign(real(r_dc))./sqrt(2);

        %Calculate BER
        BER_SC(j, i) = sum(data_sc(:) ~= data(:)) / N;
        BER_MRC(j, i) = sum(data_mrc(:) ~= data(:)) / N;
        BER_EGC(j, i) = sum(data_egc(:) ~= data(:)) / N;
        BER_DC(j, i) = sum(data_dc(:) ~= data(:)) / N;
    
    end
    
end

%Plot BER versus Es/N0 curve (Rician)
figure
semilogy(Es_N0, BER_SC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_SC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_SC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_SC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Selective Combining (Rician)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_MRC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_MRC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_MRC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_MRC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Maximal Ratio Combining (Rician)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_EGC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_EGC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_EGC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_EGC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Equal Gain Combining (Rician)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

figure
semilogy(Es_N0, BER_DC(1,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_DC(2,:),'-o','linewidth',1);
semilogy(Es_N0, BER_DC(3,:),'-o','linewidth',1);
semilogy(Es_N0, BER_DC(4,:),'-o','linewidth',1);
legend('L=1','L=2','L=3','L=4');
title('BER for QPSK with Direct Combining (Rician)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');

%Comparison of Rayleigh and Rician fading with Maximal Ratio Combining
figure
semilogy(Es_N0, BER_MRC_Rayleigh(4,:),'-o','linewidth',1);
hold on;
semilogy(Es_N0, BER_MRC(4,:),'-o','linewidth',1);
legend('Rayleigh','Rician');
title('BER for QPSK with MRC (Rayleigh & Rician)');
xlabel('Es/N0(dB)');
ylabel('Bit Error Rate(%)');