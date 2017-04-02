function [erro] = Otimizador(miks_kin, T, X, S, P, t, Pmax, YSX, YPX, esc)
	
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
	
	X0 = X(1);
	S0 = S(1);
	P0 = 0;

	%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
	Y = lsode("odefcn", [X0; S0; P0], (T = linspace(T(1), T(end), 1201 )));
	
	X1 = Y(:,1);
	S1 = Y(:,2);
	P1 = Y(:,3);
	
	X1 = X1(1:200:1201);
	S1 = S1(1:200:1201);
	P1 = P1(1:200:1201);
	
	erro = norm(X - X1) + norm(S - S1) + norm(P - P1);
	disp(erro);
end