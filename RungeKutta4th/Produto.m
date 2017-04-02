function [dPdt] = Produto(t, X, S, P, MAX_mu_X, KS, Ki, n, Pmax, YPX, esc)
    
	if esc == 1
		% Monod
		dPdt = ((MAX_mu_X * S * X)/(KS + S)) * (YPX);
	end
	if esc == 2
		% Andrews
		dPdt = (X * MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (YPX);
	end
	if esc == 3
		% Levenspiel
		dPdt = (X * (MAX_mu_X * (S/(KS + S))) * (1 - (P/Pmax)).^n) * (YPX);
	end
	if esc == 4
		% Ghose and Tyel
		dPdt = (X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax))) * (YPX);
	end
	if esc == 5
		% Jin
		dPdt = X * (MAX_mu_X * exp(-0.06035*P -0.0055*S)*(S/(KS + S))) * (YPX);
	end
end