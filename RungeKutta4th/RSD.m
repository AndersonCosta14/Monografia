function [dp] = RSD(y, yi, n)
    sum = 0;
    for i = 1:n
        sum = sum + (y(i) - yi(i))^2;
    end
    dp = sum/n;
end