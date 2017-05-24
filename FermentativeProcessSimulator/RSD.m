function [rsd] = RSD(y, yi, n)
    sum = 0;
	med = 0;
    for i = 1:n
        sum = sum + (y(i) - yi(i))^2;
		med = med + y(i);
    end
	med = med/n;
    aux = sum/n;
	
	rsd = (sqrt(aux)/med)*100;
end