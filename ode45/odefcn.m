%%
% Sistema de Equações Diferenciais
% y(1) = X(t) 
% dy(1) = dX/dt 
% y(2) = S(t) 
% dy(2) = dS/dt 
% y(3) = P(t) 
% dy(3) = dP/dt
%%

function dy = odefcn(t, y)
	
	global MAX_mu_X;
    global KS;
	global Ki;
	global n;
	global Pmax;
	global YSX;
	global YPX;
	global esc;
	
	dy = zeros(3,1);
	if esc == 1
		dy(1) = (MAX_mu_X * y(2) * y(1))/(KS + y(2));
		dy(2) = -((MAX_mu_X * y(2) * y(1))/(KS + y(2))) * (YSX);
		dy(3) = ((MAX_mu_X * y(2) * y(1))/(KS + y(2))) * (YPX);
	end
	if esc == 2
		dy(1) = y(1) * MAX_mu_X * (y(2)/(KS + y(2) + ((y(2).^2)/Ki)));
		dy(2) = -(y(1) * MAX_mu_X * (y(2)/(KS + y(2) + ((y(2).^2)/Ki)))) * (YSX);
		dy(3) = (y(1) *MAX_mu_X *(y(2)/(KS + y(2) + ((y(2).^2)/Ki)))) * (YPX);
	end
	if esc == 3
		dy(1) = y(1) * (MAX_mu_X * (y(2)/(KS + y(2)))) * (1 - (y(3)/Pmax)).^n;
		dy(2) = -(y(1) * (MAX_mu_X* (y(2)/(KS + y(2)))) * (1 - (y(3)/Pmax)).^n) * (YSX);
		dy(3) = (y(1) * (MAX_mu_X * (y(2)/(KS + y(2)))) * (1 - (y(3)/Pmax)).^n) * (YPX);
	end
	if esc == 4
		dy(1) = y(1) * (MAX_mu_X * (y(2)/(KS + y(2) + ((y(2).^2)/Ki)))) * (1 - (y(3)/Pmax));
		dy(2) = -(y(1) * (MAX_mu_X * (y(2)/(KS + y(2) + ((y(2).^2)/Ki))))) * (1 - (y(3)/Pmax)) * (YSX);
		dy(3) = (y(1) * (MAX_mu_X * (y(2)/(KS + y(2) + ((y(2).^2)/Ki)))) * (1 - (y(3)/Pmax))) * (YPX);
	end
	if esc == 5
		dy(1) = y(1) * MAX_mu_X * exp(-0.06035*y(3) -0.0055*y(2)) * (y(2)/(KS + y(2)));
		dy(2) = -y(1) * (MAX_mu_X * exp(-0.06035*y(3) -0.0055*y(2)) * (y(2)/(KS + y(2)))) * (YSX);
		dy(3) = y(1) * (MAX_mu_X * exp(-0.06035*y(3) -0.0055*y(2))*(y(2)/(KS + y(2)))) * (YPX);
	end
	
	options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
end