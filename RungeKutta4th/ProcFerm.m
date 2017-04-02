%% MONOGRAFIA: Código para Simulação do Processo Fermentativo do Etanol
%  ALUNO:      Anderson da Silva Costa
%  MATRÍCULA:  358551
%  ORIENTADOR: Michael Ferreira de Souza
%%
function ProcFerm()
    run = 1;
    while (run == 1)
        clc;
        clear all;
        choice = menu('Concentração Inicial de Substrato','70 g/L','80 g/L','90 g/L','110 g/L','130 g/L','170 g/L');
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
		
		X0 = X(1);
		S0 = S(1);
		P0 = 0;
				
		% Fatores de Conversão
		YSX = (S(1) - S(end))/(X(end) - X(1));
		YPX = (P(end) - P(1))/(X(end) - X(1));
				
		% Interpolação de Splines
        polyX = @(t) interp1(T, X, t, 'spline');
        polyS = @(t) interp1(T, S, t, 'spline');
        polyP = @(t) interp1(T, P, t, 'spline');
                
		% Valores máximos de Concentração
        [~, Xmax] = fminbnd(@(t) -polyX(t), T(1), T(end));
        Xmax = -Xmax;
               
        [~, Smax] = fminbnd(@(t) -polyS(t), T(1), T(end));
        Smax = -Smax;
                
        [~, Pmax] = fminbnd(@(t) -polyP(t), T(1), T(end));
        Pmax = -Pmax;
                
		%% Calculo dos Parâmetros dos Modelos
				
		% Discretização do tempo descartando os tempos cuja concentração de Substrato é nula
		h     = 0.001;
		t     = T(1):h:T(min(find(S,1,'last') + 1, length(T)));
				
		% Valores de concentração de Biomassa no Tempo discretizado
		X_mu  = polyX(t);
				
		% Aproximação da Derivada por diferenças finitas
		dX_mu = diff(X_mu)/h;
				
		% Calculo da Velocidade Específica de Transformação da Biomassa (mu_X)
		mu_X  = dX_mu./X_mu(1:end-1);
		t     = t(1:end-1);
		X_mu  = X_mu(1:end-1);
				
		% Valores da concentração de Substrato no Tempo
		S_mu  = polyS(t);
				
		% Determinação do valor maximo de mu_X
		[MAX_mu_X, ~] = max(mu_X);
		aux = MAX_mu_X/2; % Metade do Valor Máximo de mu_X
				
		% Determinar a posição de aux (Valor Máximo dividido por 2)
		for i=1:length(mu_X)
			if abs(mu_X(i) - aux) <= 1e-4
				pos = i;
			end
		end
				
		% Constante de Saturação
		KS = S_mu(pos);
				
		% Constante de Inibição pelo Produto e parâmetro n
		Ki = 1;
		n  = 1;
				
		modelo = menu('Modelo Matemático', 'Monod', 'Andrews', 'Levenspiel', 'Ghose and Tyel', 'Jin');
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
        % Otimização dos Parâmetros				
		coef = lsqnonlin(@(miks_kin)calc_coef(miks_kin, mu_X, X, polySt2, polyPt2, Pmax, esc),[MAX_mu_X, KS, Ki, n]);
				
        x = lsqnonlin(@(miks)Otimizador(miks, T, X, S, P, t, Pmax, YSX, YPX, esc), coef);
						
		mu_MAX = x(1);
		KS     = x(2);
		Ki     = x(3);
		n      = x(4);
						
		% Resolução do Sistema de EDO's com os melhores parâmetros possíveis.
		[t1, X1, S1, P1] = RungeKutta4th('Biomassa', 'Substrato', 'Produto', T(1), T(end), 0.001, X0, S0, P0, mu_MAX, KS, Ki, n, Pmax, YSX, YPX, esc);
		
		time = cputime-t;
		
		X1 = X1';
		S1 = S1';
		P1 = P1';
		auxX = X1;
		auxS = S1;
		auxP = P1;
		X1 = X1(1:2000:T(end)*1000+1);
		S1 = S1(1:2000:T(end)*1000+1);
		P1 = P1(1:2000:T(end)*1000+1);
		
		% Cálculo do Desvio Padrão Residual
		dpX = RSD(X, X1, length(X));
		dpS = RSD(S, S1, length(S));
		dpP = RSD(P, P1, length(P));
		
		% Teste de Confiança utilizando nível de Randomicidade
		Z1 = Randomic(X, X1);
		Z2 = Randomic(S, S1);
		Z3 = Randomic(P, P1);
		
		% Resultados
		fileID = fopen('Resultados.txt','wt');
		switch esc
			case 1 % Monod
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Monod com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 2 % Andrews
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Andrews com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Constante de Inibição:                ', Ki);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 3 % Levenspiel
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Levenspiel com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Expoente n:                           ', n);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 4 % Ghoose and Tyagi
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Ghose and Tyagi com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Constante de Inibição:                ', Ki);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 5 % Jin
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Jin com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, 'g/L');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
		end
		fclose(fileID);
		
		% Plotagem dos Dados
		fig = figure; set(fig,'Position',[320 640 640 320]);
				
		polyX = @(w) interp1(t1, auxX, w, 'spline');
		polyXt = polyX(T);
		subplot(1,3,1); hold on;  box on;
		plot(T, polyXt, '-r', 'LineWidth', 2);
		plot(T, X, '*');
		title('Biomassa');
		ylabel('Concentracao de Biomassa - X(t)');
		xlabel('Tempo - t');

		polyS = @(w) interp1(t1, auxS, w, 'spline');
		polySt = polyS(T);
		subplot(1,3,2); hold on;  box on;
		plot(T, polySt, '-r', 'LineWidth', 2);
		plot(T, S, '*');
		title('Substrato');
		ylabel('Concentracao de Substrato - S(t)');
		xlabel('Tempo - t');

		polyP = @(w) interp1(t1, auxP, w, 'spline');
		polyPt = polyP(T);
		subplot(1,3,3); hold on;  box on;
		plot(T, polyPt, '-r', 'LineWidth', 2);
		plot(T, P, '*');
		title('Produto');
		ylabel('Concentracao de Produto - P(t)');
		xlabel('Tempo - t');
					
		drawnow
		
        stop = menu('Deseja realizar um novo teste?','Sim','Não');
        if stop == 1
            run = 1;
        end
        if stop == 2
            run = 0;
        end
    end
end