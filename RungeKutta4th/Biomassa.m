function [dXdt] = Biomassa(t, X, S, P, MAX_mu_X, KS, Ki, n, Pmax, esc)
	
	if esc == 1
		% Monod
		dXdt = (MAX_mu_X * S * X)/(KS + S);
	end
	if esc == 2
		% Andrews
		dXdt = X * MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)));
	end
	if esc == 3
		% Levenspiel
		dXdt = X * (MAX_mu_X * (S/(KS + S))) * (1 - (P/Pmax)).^n;
	end
	if esc == 4
		% Ghose and Tyel
		dXdt = X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax));
	end
	if esc == 5
		% Jin
		dXdt = X * MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S));
	end
end