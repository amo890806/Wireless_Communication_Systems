function xn = Newton_Method(yn, x0)

n = 200;
error = 1e-6;

x = x0;

for i=1:n
    [y, y_diff] = Erlang_B(x, x0);
    y_error = y - yn;
    if(abs(y_error) < error)
        break
    end
    x = x - y_error/y_diff;
    
    if x < 0
        x = x0*randn(1);
    end
end

xn = x;

end