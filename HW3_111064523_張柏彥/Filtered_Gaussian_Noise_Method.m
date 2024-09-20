close all;
clear;
clc;

num = 3;
fmT_num = [0.01, 0.1, 0.5];
sample_point = 50000;
t_over_T = 300;
Omega = 1;


for n=1:num
    
    fmT = fmT_num(n);
    fm = fmT;
    tau = 10 / fm;
    
    %Set 3dB point at fm/4, eta is the coeffient of the first-order LPF
    tmp = 2 - cos(pi*fmT/2);
    eta = tmp - sqrt(tmp.^2-1);

    %Variance of Zero Mean Gaussian noise source
    var = ((1+eta)/(1-eta)) * (Omega/2);
    sigma = sqrt(var);
    pd = makedist('Normal','mu',0,'sigma',sigma);
    w1 = random(pd,sample_point,1);
    w2 = random(pd,sample_point,1);

    %Generate gI, gQ, r
    gI_init = 1/sqrt(2);
    gQ_init = 1/sqrt(2);
    gI = zeros(sample_point, 1);
    gQ = zeros(sample_point, 1);
    r = zeros(sample_point, 1);
    for i=1:sample_point
        if i==1
            gI(i,:) = eta*gI_init + (1-eta)*w1(i,:);
            gQ(i,:) = eta*gQ_init + (1-eta)*w2(i,:);
        else
            gI(i,:) = eta*gI(i-1,:) + (1-eta)*w1(i,:);
            gQ(i,:) = eta*gQ(i-1,:) + (1-eta)*w2(i,:);
        end
        r(i,:) = gI(i,:) + 1j*gQ(i,:);
    end

    %Calculate Envelope
    envelope = zeros(t_over_T, 1);
    for i=1:t_over_T
        envelope(i,:) = sqrt(gI(i,:).^2 + gQ(i,:).^2);
    end

    %Calculate Autocorrelation
    phi = zeros(tau+1, 1);
    r1 = r;
    r2 = r;
    phi_tmp = r1 .* conj(r2);
    phi(1,:) = sum(phi_tmp / (norm(r1)*norm(r2)));
    for i=1:tau
        r1 = [r1; 0];
        r2 = [0; r2];
        phi_tmp = r1 .* conj(r2);
        phi(i+1,:) = sum(phi_tmp / (norm(r1)*norm(r2)));
    end
    
    x_envelope = [1:1:t_over_T]';
    if n==1
        x1_autocorrelation = linspace(0, fm*tau, tau+1)';
        y1_envelope = 10*log10(envelope);
        y1_autocorrelation = real(phi);
    elseif n==2
        x2_autocorrelation = linspace(0, fm*tau, tau+1)';
        y2_envelope = 10*log10(envelope);
        y2_autocorrelation = real(phi);
    else
        x3_autocorrelation = linspace(0, fm*tau, tau+1)';
        y3_envelope = 10*log10(envelope);
        y3_autocorrelation = real(phi);
    end
    
    
end


%Plot Envelope
figure(1)
hold on
plot(x_envelope, y1_envelope, '-b')
plot(x_envelope, y2_envelope, '-g')
plot(x_envelope, y3_envelope, '-r')
hold off
ylim([-14 4]);
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('Filtered Gaussian method Envelope Level');
xlabel('Time, t/T');
ylabel('Envelope Level');

%Plot AutoCorrelation
figure(2)
hold on
plot(x1_autocorrelation, y1_autocorrelation, '-b')
plot(x2_autocorrelation, y2_autocorrelation, '-g')
plot(x3_autocorrelation, y3_autocorrelation, '-r')
hold off
ylim([-0.2 1]);
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('Filtered Gaussian method Autocorrelation');
xlabel('Time Delay, fm*tau');
ylabel('Autocorrelation');


