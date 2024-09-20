close all;
clear;
clc;

num = 3;
fmT_num = [0.01, 0.1, 0.5];
sample_point = 50000;
t_over_T = 300;
Omega = 1;
M = 8;
N = 2*(2*M+1);
n = [1:1:N];
theta_n = 2*pi*n / N;
theta_n = theta_n(1:M);
beta_n = pi*n/M;
beta_n = beta_n(1:M);
alpha = 0;

for n=1:num

    fmT = fmT_num(n);
    fm = fmT;
    fn = fm * cos(theta_n);
    tau = 10 / fm;

    cos_beta_n = cos(beta_n);
    sin_beta_n = sin(beta_n);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);

    gI = zeros(sample_point, 1);
    gQ = zeros(sample_point, 1);
    r = zeros(sample_point, 1);
    for t=0:sample_point
        cos_fnt = cos(2*pi*fn*t);
        cos_fmt = cos(2*pi*fm*t);
        gI(t+1,:) = sqrt(2)*(2*sum(cos_beta_n.*cos_fnt) + sqrt(2)*cos_alpha*cos_fmt);
        gQ(t+1,:) = sqrt(2)*(2*sum(sin_beta_n.*cos_fnt) + sqrt(2)*sin_alpha*cos_fmt);
        r(t+1,:) = gI(t+1,:) + 1j*gQ(t+1,:);
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
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('M=8, Jake''s method Envelope Level');
xlabel('Time, t/T');
ylabel('Envelope Level');

%Plot AutoCorrelation
figure(2)
hold on
plot(x1_autocorrelation, y1_autocorrelation, '-b')
plot(x2_autocorrelation, y2_autocorrelation, '-g')
plot(x3_autocorrelation, y3_autocorrelation, '-r')
hold off
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('M=8, Jake''s method Autocorrelation');
xlabel('Time Delay, fm*tau');
ylabel('Autocorrelation');

num = 3;
fmT_num = [0.01, 0.1, 0.5];
sample_point = 50000;
t_over_T = 300;
Omega = 1;
M = 16;
N = 2*(2*M+1);
n = [1:1:N];
theta_n = 2*pi*n / N;
theta_n = theta_n(1:M);
beta_n = pi*n/M;
beta_n = beta_n(1:M);
alpha = 0;

for n=1:num

    fmT = fmT_num(n);
    fm = fmT;
    fn = fm * cos(theta_n);
    tau = 10 / fm;

    cos_beta_n = cos(beta_n);
    sin_beta_n = sin(beta_n);
    cos_alpha = cos(alpha);
    sin_alpha = sin(alpha);

    gI = zeros(sample_point, 1);
    gQ = zeros(sample_point, 1);
    r = zeros(sample_point, 1);
    for t=0:sample_point
        cos_fnt = cos(2*pi*fn*t);
        cos_fmt = cos(2*pi*fm*t);
        gI(t+1,:) = sqrt(2)*(2*sum(cos_beta_n.*cos_fnt) + sqrt(2)*cos_alpha*cos_fmt);
        gQ(t+1,:) = sqrt(2)*(2*sum(sin_beta_n.*cos_fnt) + sqrt(2)*sin_alpha*cos_fmt);
        r(t+1,:) = gI(t+1,:) + 1j*gQ(t+1,:);
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
figure(3)
hold on
plot(x_envelope, y1_envelope, '-b')
plot(x_envelope, y2_envelope, '-g')
plot(x_envelope, y3_envelope, '-r')
hold off
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('M=16, Jake''s method Envelope Level');
xlabel('Time, t/T');
ylabel('Envelope Level');

%Plot AutoCorrelation
figure(4)
hold on
plot(x1_autocorrelation, y1_autocorrelation, '-b')
plot(x2_autocorrelation, y2_autocorrelation, '-g')
plot(x3_autocorrelation, y3_autocorrelation, '-r')
hold off
ylim([-0.6 1]);
legend('fmT=0.01','fmT=0.1','fmT=0.5');
title('M=16, Jake''s method Autocorrelation');
xlabel('Time Delay, fm*tau');
ylabel('Autocorrelation');

