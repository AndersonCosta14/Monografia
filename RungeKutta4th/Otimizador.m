function [erro] = Otimizador(miks_kin, T, X, S, P, t, Pmax, YSX, YPX, esc)
	MAX_mu_X = miks_kin(1);
	KS       = miks_kin(2);
	Ki       = miks_kin(3);
	n        = miks_kin(4);
	
	X0 = X(1);
	S0 = S(1);
	P0 = 0;

	[t, X1, S1, P1] = RungeKutta4th('Biomassa', 'Substrato', 'Produto', T(1), T(end), 0.001, X0, S0, P0, MAX_mu_X, KS, Ki, n, Pmax, YSX, YPX, esc);
	
	X1 = X1';
	S1 = S1';
	P1 = P1';
	X1 = X1(1:2000:T(end)*1000+1);
	S1 = S1(1:2000:T(end)*1000+1);
	P1 = P1(1:2000:T(end)*1000+1);
	
	erro = norm(X - X1) + norm(S - S1) + norm(P - P1);
	disp(erro);
end