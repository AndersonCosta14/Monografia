function [dSdt] = Substrato(t, X, S, P, MAX_mu_X, KS, Ki, n, Pmax, YSX, esc)

	if esc == 1
		% Monod
		dSdt = -((MAX_mu_X * S * X)/(KS + S)) * (YSX);
	end
	if esc == 2
		% Andrews
		dSdt = -(X * MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (YSX);
	end
	if esc == 3
		% Levenspiel
		dSdt = -(X * (MAX_mu_X * (S/(KS + S))) * (1 - (P/Pmax)).^n) * (YSX);
	end
	if esc == 4
		% Ghose and Tyel
		dSdt = -(X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki))))) * (1 - (P/Pmax)) * (YSX);
	end
	if esc == 5
		% Jin
		dSdt = -X * (MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S))) * (YSX);
	end
end