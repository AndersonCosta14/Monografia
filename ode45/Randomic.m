function [Z] = Randomic(y, yi)
    %y : dados experimentais
    %yi: dados do modelo
    
    N1 = 0;
    N2 = 0;
    Res = zeros(1, length(y));
    R = 0;
    
    for i = 1:length(y)
        Res(i) = yi(i) - y(i); 
        if(Res(i) > 0)
            N1 = N1 + 1;
        end
        if(Res(i) < 0)
            N2 = N2 + 1;
        end
    end
    for i = 1:length(y)-1
        if(Res(i)*Res(i+1) < 0)
            R = R + 1;
        end
    end
    mediaR = ((2*N1*N2)/(N1 + N2)) + 1;
    sigmaR = sqrt((2*N1*N2*(2*N1*N2 - N1 - N2))/(((N1 + N2)^2 )*(N1 + N2 - 1)));
    
    Z = (R - mediaR)/sigmaR;
end