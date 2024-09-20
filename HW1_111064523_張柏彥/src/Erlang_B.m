function [B, B_diff] = Erlang_B(rho, m)

sum = 0;
sum_u_diff = 0;
sum_v_diff = 0;
for k = 0:m

    prod = 1;
    for i = k+1 : m
        prod = prod * i / rho;
    end

    temp = prod;
    temp_diff = (k-m) * prod / rho;

    sum = sum + temp;
    sum_u_diff = sum_u_diff + temp_diff;
    sum_v_diff = sum_v_diff + temp;

end

B = 1 / sum;
B_diff = -sum_u_diff / sum_v_diff^2;
    
end