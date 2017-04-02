%% MONOGRAFIA: Cdigo para Simulao do Processo Fermentativo do Etanol
%  ALUNO:      Anderson da Silva Costa
%  MATRCULA:  358551
%  ORIENTADOR: Michael Ferreira de Souza
%%
function ProcFerm()
    run = 1;
    while (run == 1)
        clc;
        clear all;
        close all;
        choice = menu('Concentrao Inicial de Substrato','70 g/L','80 g/L','90 g/L','110 g/L','130 g/L','170 g/L');
        T = 0:2:12;
        T = T';
		
        switch choice
            case 1
                [X, S, P] = DadosExp(70)
				Inicial = 70;
            case 2
                [X, S, P] = DadosExp(80)
				Inicial = 80;
            case 3
                [X, S, P] = DadosExp(90)
				Inicial = 90;
            case 4
                [X, S, P] = DadosExp(110)
				Inicial = 110;
            case 5
                [X, S, P] = DadosExp(130)
				Inicial = 130;
            case 6
                [X, S, P] = DadosExp(170)
				Inicial = 170;
        end
		
		global MAX_mu_X;
        global KS;
		global Ki;
		global n;
		global Pmax;
		global YSX;
		global YPX;
		global esc;
		
		% Condies Iniciais
		X0 = X(1);
		S0 = S(1);
		P0 = 0;
				
		% Fatores de Converso
		YSX = (S(1) - S(end))/(X(end) - X(1));
		YPX = (P(end) - P(1))/(X(end) - X(1));
				
		% Interpolao de Splines
        polyX = @(t) interp1(T, X, t, 'spline');
        polyS = @(t) interp1(T, S, t, 'spline');
        polyP = @(t) interp1(T, P, t, 'spline');
                
		% Valores mximos de Concentrao
        [~, Xmax] = fminbnd(@(t) -polyX(t), T(1), T(end));
        Xmax = -Xmax;
               
        [~, Smax] = fminbnd(@(t) -polyS(t), T(1), T(end));
        Smax = -Smax;
                
        [~, Pmax] = fminbnd(@(t) -polyP(t), T(1), T(end));
        Pmax = -Pmax;
                
		%% Calculo dos Parmetros dos Modelos
				
		% Discretizao do tempo descartando os tempos cuja concentrao de Substrato  nula
		h     = 0.001;
		t     = T(1):h:T(min(find(S,1,'last') + 1, length(T)));
				
		% Valores de concentrao de Biomassa no Tempo discretizado
		X_mu  = polyX(t);
				
		% Aproximao da Derivada por diferenas finitas
		dX_mu = diff(X_mu)/h;
				
		% Calculo da Velocidade Especfica de Transformao da Biomassa (mu_X)
		mu_X  = dX_mu./X_mu(1:end-1);
		t     = t(1:end-1);
		X_mu  = X_mu(1:end-1);
				
		% Valores da concentrao de Substrato no Tempo
		S_mu  = polyS(t);
				
		% Determinao do valor maximo de mu_X
		[MAX_mu_X, ~] = max(mu_X);
		aux = MAX_mu_X/2; % Metade do Valor Mximo de mu_X
				
		% Determinar a posio de aux (Valor Mximo dividido por 2)
		for i=1:length(mu_X)
			if abs(mu_X(i) - aux) <= 1e-4
				pos = i;
			end
		end
				
		% Constante de Saturao
		KS = S_mu(pos);
				
		% Constante de Inibio pelo Produto e parmetro n
		Ki = 1;
		n  = 1;
				
		modelo = menu('Modelo Matemtico', 'Monod', 'Andrews', 'Levenspiel', 'Ghose and Tyel', 'Jin');
		polySt2   = polyS(t);
		polyPt2   = polyP(t);
		switch modelo
			case 1
				esc = 1;
			case 2
				esc = 2;
			case 3
				esc = 3;
			case 4
				esc = 4;
			case 5
				esc = 5;
        end
		
		t = cputime;
        % Otimizao dos Parmetros	
		
		coef = fminsearch(@(miks_kin)calc_coef(miks_kin, mu_X, X, polySt2, polyPt2, Pmax, esc),[MAX_mu_X, KS, Ki, n]);
				
        x = fminsearch(@(miks)Otimizador(miks, T, X, S, P, t, Pmax, YSX, YPX, esc), coef);
						
		mu_MAX = x(1);
		KS     = x(2);
		Ki     = x(3);
		n      = x(4);
		
		% Resoluo do Sistema de EDO's com os melhores parmetros possveis.

		%options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
		T1 = linspace(T(1), T(end), 1201);
		Y = lsode("odefcn", [X0; S0; P0], T1);
			
		X1 = Y(:,1);
		S1 = Y(:,2);
		P1 = Y(:,3);
				
		time = cputime-t;
		
		Xaux = X1(1:200:1201);
	          	Saux = S1(1:200:1201);
		Paux = P1(1:200:1201);
		
		% Clculo do Desvio Padro Residual
		dpX = RSD(X, Xaux, length(X));
		dpS = RSD(S, Saux, length(S));
		dpP = RSD(P, Paux, length(P));
		
		% Teste de Confiana utilizando nvel de Randomicidade
		Z1 = Randomic(X, Xaux);
		Z2 = Randomic(S, Saux);
		Z3 = Randomic(P, Paux);
		
		% Resultados
		fileID = fopen('Resultados.txt','wt');
		switch esc
			case 1 % Monod
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Monod com concentrao inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu mximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturao:               ', KS);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulao:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Produto:    ', Z3);
			case 2 % Andrews
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Andrews com concentrao inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu mximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturao:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Constante de Inibio:                ', Ki);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulao:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Produto:    ', Z3);
			case 3 % Levenspiel
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Levenspiel com concentrao inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu mximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturao:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Expoente n:                           ', n);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulao:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Produto:    ', Z3);
			case 4 % Ghoose and Tyagi
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Ghose and Tyagi com concentrao inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu mximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturao:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Constante de Inibio:                ', Ki);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulao:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Produto:    ', Z3);
			case 5 % Jin
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Jin com concentrao inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentrao Mxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu mximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturao:               ', KS);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulao:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padro Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nvel de Randomicidade do Produto:    ', Z3);
		end
		fclose(fileID);

		% Plotagem dos Dados
		clf ();
		fig = figure(); 
		
		subplot(1,3,1); 
		hold on;  
		plot(T1, X1, '-r', 'LineWidth', 2);
		plot(T,X,'d');
		xlabel('Tempo (h)');
		ylabel('Concentracao de Biomassa (g/L)');
		
		subplot(1,3,2);
		hold on;  
		plot(T1, S1, '-r', 'LineWidth', 2);
		plot(T,S,'d');
		xlabel('Tempo (h)');
		ylabel('Concentracao de Substrato (g/L)');
		
		subplot(1,3,3); 
		hold on;  
		plot(T1, P1, '-r', 'LineWidth', 2);
		plot(T,P,'d');
		xlabel('Tempo (h)');
		ylabel('Concentracao de Produto (g/L)');
		
		legend('Numerica','Empirica','Location','southeast')
		print -djpg fig.jpg
		
        stop = menu('Deseja realizar um novo teste?','Sim','Nao');
        if stop == 1
            run = 1;
        end
        if stop == 2
            run = 0;
        end
    end
end