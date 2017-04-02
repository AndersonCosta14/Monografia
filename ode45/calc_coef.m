function [erro] = calc_coef(miks_kin, mu_X, X, S, P, Pmax, esc)
	
	global MAX_mu_X;
    global KS;
	global Ki;
	global n;
	global Pmax;
	global YSX;
	global YPX;
	global esc;
	
	MAX_mu_X = miks_kin(1);
	KS       = miks_kin(2);
	Ki       = miks_kin(3);
	n        = miks_kin(4);
	
	if esc == 1
		% Monod
		resp = (MAX_mu_X * S)/(KS + S);
	end
	if esc == 2
		% Andrews
		resp = MAX_mu_X*(S/(KS + S + ((S.^2)/Ki)));
	end
	if esc == 3
		% Levenspiel
		resp = (MAX_mu_X*(S/(KS + S))) * (1 - (P/Pmax)).^n;
	end
    if esc == 4
		% Ghose and Tyel
		resp = (MAX_mu_X*(S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax));
	end
	if esc == 5
		%Jin
		resp = MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S));
	end
	
	erro = norm(resp - mu_X);
end