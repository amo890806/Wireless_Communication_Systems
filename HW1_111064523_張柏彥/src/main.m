close all;
clear;
clc;

ans1_20 = zeros(20, 4);
ans200_220 = zeros(21, 4);
Block_Rate = [0.01 0.03 0.05 0.1];
for i = 1:4
    for m = 1:20
        rho = Newton_Method(Block_Rate(i), m);
        ans1_20(m,i) = rho;
    end
end

for i = 1:4
    n = 1;
    for m = 200:220
        rho = Newton_Method(Block_Rate(i), m);
        ans200_220(n,i) = rho;
        n = n + 1;
    end
end