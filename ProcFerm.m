%% 
% UNIVERSIDADE FEDERAL DO CEARÁ
% DEPARTAMENTO DE ESTATÍSTICA E MATEMÁTICA APLICADA
% Monografia: Código para Simulação do Processo Fermentativo do Etanol
% Aluno: Anderson da Silva Costa
% Matrícula: 358551
% Orientador: Michael Ferreira de Souza
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
		
		global MAX_mu_X;
        global KS;
		global Ki;
		global n;
		global Pmax;
		global YSX;
		global YPX;
		global esc;
		
		% Condições Iniciais
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
                
		%% Cálculo dos Parâmetros dos Modelos
				
		% Discretização do tempo descartando os tempos cuja concentração de Substrato é nula
		h     = 0.001;
		t     = T(1):h:T(min(find(S,1,'last') + 1, length(T)));
				
		% Valores de concentração de Biomassa no Tempo discretizado
		X_mu  = polyX(t);
				
		% Aproximação da Derivada por diferenÃ§as finitas
		dX_mu = diff(X_mu)/h;
				
		% Cálculo da Velocidade Específica de Transformação da Biomassa (mu_X)
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
				
		modelo = menu('Modelo Matemático', 'Monod', 'Andrews', 'Levenspiel', 'Ghose e Thyagi', 'Jin, Chiang e Wang');
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
		
		% Otimização dos Parâmetros	
		coef = lsqnonlin(@(miks_kin)calc_coef(miks_kin, mu_X, X, polySt2, polyPt2, Pmax, esc),[MAX_mu_X, KS, Ki, n]);
		
		metodo = menu('Métodos Numéricos', 'ode45', 'Runge-Kutta 4th');
		switch metodo
			case 1
			    t = cputime;
				
				x = lsqnonlin(@(miks)OtimizadorODE45(miks, T, X, S, P, t, Pmax, YSX, YPX, esc), coef);
						
				mu_MAX = x(1);
				KS     = x(2);
				Ki     = x(3);
				n      = x(4);
		
				% Resolução do Sistema de EDO's utilizando a ode45 com os melhores parâmetros possíveis.

				options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
				[T1, Y] = ode45(@odefcn, [T(1):2:T(end)], [X0 S0 P0], options);
			
				X1 = Y(:,1);
				S1 = Y(:,2);
				P1 = Y(:,3);
				
				time = cputime-t;
			case 2
				t = cputime;
				x = lsqnonlin(@(miks)OtimizadorRK4th(miks, T, X, S, P, t, Pmax, YSX, YPX, esc), coef);
						
				mu_MAX = x(1);
				KS     = x(2);
				Ki     = x(3);
				n      = x(4);
						
				% Resolução do Sistema de EDO's através o Runge-Kutta 4th com os melhores parâmetros possíveis.
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
		end
		
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
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, '%');
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
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, '%');
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
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, '%');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 4 % Ghoose e Thyagi
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Ghose e Thyagi com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f\n', 'Constante de Inibição:                ', Ki);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, '%');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
			case 5 % Jin, Chiang e Wang
				fprintf(fileID, '%12s %d %3s\n\n', 'Modelo de Jin, Chiang e Wang com concentração inicial de', Inicial, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Biomassa:      ', Xmax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Substrato:     ', Smax, 'g/L');
				fprintf(fileID, '%12s %f %3s\n', 'Concentração Máxima de Produto:       ', Pmax, 'g/L');
				fprintf(fileID, '%12s %f\n', 'Mu máximo:                            ', mu_MAX);
				fprintf(fileID, '%12s %f\n', 'Constante de Saturação:               ', KS);
				fprintf(fileID, '%12s %f %6s\n', 'Tempo de Simulação:                   ', time, 'seconds');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual da Biomassa:   ', dpX, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Substrato:  ', dpS, '%');
				fprintf(fileID, '%12s %f %6s\n', 'Desvio Padrão Residual do Produto:    ', dpP, '%');
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade da Biomassa:   ', Z1);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Substrato:  ', Z2);
				fprintf(fileID, '%12s %f\n', 'Nível de Randomicidade do Produto:    ', Z3);
		end
		fclose(fileID);
		
		% Plotagem dos Dados
		fig = figure(); 
		set(fig,'Position',[320 640 640 320]);
		
		subplot(1,3,1); 
        hold on;  
		plot(T, X1, '--b', 'LineWidth', 1);
		plot(T, X, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'b');
		title('Biomassa');
		ylabel('Concentração de Biomassa (g/L)');
		xlabel('Tempo (h)');
		legend('X (Numérico)', 'X (Empírico)', 'Location', 'southeast')

		subplot(1,3,2); 
        hold on;  
		plot(T, S1, '--k', 'LineWidth', 1);
		plot(T, S, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
		title('Substrato');
		ylabel('Concentração de Substrato (g/L)');
		xlabel('Tempo (h)');
		legend('S (Numérico)', 'S (Empírico)', 'Location', 'northeast')
		
		subplot(1,3,3); 
        hold on;  
		plot(T, P1, '--r', 'LineWidth', 1);
		plot(T, P, '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'r');
		title('Produto');
		ylabel('Concentração de Produto (g/L)');
		xlabel('Tempo (h)');
		legend('P (Numérico)', 'P (Empírico)', 'Location', 'southeast')
		
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

function [X, S, P] = DadosExp(S0)
    switch S0
        case 70
            fprintf('Valores de concentração para Substrato Inicial de 70 g/L.');
            X = [3.680; 4.700; 5.520; 5.800; 6.220; 5.840; 5.910];
            S = [72.405; 45.985; 17.420; 0.000; 0.000; 0.000; 0.000];
            P = [0.000; 11.050; 24.375; 30.520; 30.950; 30.785; 31.205];
        case 80
            fprintf('Valores de concentração para Substrato Inicial de 80 g/L.');
            X = [3.920; 4.840; 6.210; 6.850; 7.140; 6.930; 6.910];
            S = [81.170; 56.140; 23.500; 8.830; 0.000; 0.000; 0.000];
            P = [0.000; 10.830; 21.400; 31.850; 34.260; 34.680; 34.650];
        case 90
            fprintf('Valores de concentração para Substrato Inicial de 90 g/L.');
            X = [3.680; 4.710; 5.450; 6.260; 6.390; 6.550; 6.420];
            S = [93.850; 73.660; 41.820; 16.030; 0.000; 0.000; 0.000];
            P = [0.000; 10.450; 24.290; 35.820; 43.560; 43.130; 44.130];
        case 110
            fprintf('Valores de concentração para Substrato Inicial de 110 g/L.');
            X = [4.210; 5.400; 6.060; 6.570; 6.770; 6.440; 6.910];
            S = [114.910; 86.420; 55.580; 30.300; 16.290; 0.000; 0.000];
            P = [0.000; 15.650; 30.400; 42.270; 46.900; 51.720; 53.540];
        case 130
            fprintf('Valores de concentração para Substrato Inicial de 130 g/L.');
            X = [4.170; 4.710; 5.490; 6.220; 6.690; 6.920; 7.010];
            S = [131.990; 110.440; 81.600; 50.040; 23.810; 0.000; 0.000];
            P = [0.000; 8.550; 22.240; 38.010; 49.510; 55.400; 57.880];
        case 170
            fprintf('Valores de concentração para Substrato Inicial de 170 g/L.');
            X = [4.280; 4.650; 5.420; 6.520; 7.470; 8.000; 8.220;];
            S = [170.060; 149.920; 123.590; 93.920; 66.210; 43.290; 24.700];
            P = [0.000; 5.990; 15.650; 28.550; 42.570; 54.780; 63.510];
        otherwise
            disp('Concentração Inicial inválida');
    end
end

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
		% Ghose e Thyagi
		dXdt = X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax));
	end
	if esc == 5
		% Jin, Chiang e Wang
		dXdt = X * MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S));
	end
end

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
		% Ghose e Thyagi
		dSdt = -(X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki))))) * (1 - (P/Pmax)) * (YSX);
	end
	if esc == 5
		% Jin, Chiang e Wang
		dSdt = -X * (MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S))) * (YSX);
	end
end

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
		% Ghose e Thyagi
		dPdt = (X * (MAX_mu_X * (S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax))) * (YPX);
	end
	if esc == 5
		% Jin, Chiang e Wang
		dPdt = X * (MAX_mu_X * exp(-0.06035*P -0.0055*S)*(S/(KS + S))) * (YPX);
	end
end

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

function [t, x, y, z] = RungeKutta4th(EDO1, EDO2, EDO3, a, b, h, x1, y1, z1, par1, par2, par3, par4, par5, par6, par7, par8);
	t(1) = a;
	x(1) = x1;
	y(1) = y1;
	z(1) = z1;

	n    = (b-a)/h;

	for i = 1:n
		t(i+1) = t(i) + h;
		tm     = t(i) + h/2;
		Kx1    = feval(EDO1, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par8);
		Ky1    = feval(EDO2, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par6, par8);
		Kz1    = feval(EDO3, t(i), x(i), y(i), z(i), par1, par2, par3, par4, par5, par7, par8);
		Kx2    = feval(EDO1, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par8);
		Ky2    = feval(EDO2, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par6, par8);
		Kz2    = feval(EDO3, tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2, z(i) + Kz1*h/2, par1, par2, par3, par4, par5, par7, par8);
		Kx3    = feval(EDO1, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par8);
		Ky3    = feval(EDO2, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par6, par8);
		Kz3    = feval(EDO3, tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2, z(i) + Kz2*h/2, par1, par2, par3, par4, par5, par7, par8);
		Kx4    = feval(EDO1, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par8);
		Ky4    = feval(EDO2, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par6, par8);
		Kz4    = feval(EDO3, t(i+1), x(i) + Kx3*h, y(i) + Ky3*h, z(i) + Kz3*h, par1, par2, par3, par4, par5, par7, par8);
		x(i+1) = x(i) + (Kx1 + 2*Kx2 + 2*Kx3 + Kx4)*h/6;
		y(i+1) = y(i) + (Ky1 + 2*Ky2 + 2*Ky3 + Ky4)*h/6;
		z(i+1) = z(i) + (Kz1 + 2*Kz2 + 2*Kz3 + Kz4)*h/6;
	end

end

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
		% Ghose e Thyagi
		resp = (MAX_mu_X*(S/(KS + S + ((S.^2)/Ki)))) * (1 - (P/Pmax));
	end
	if esc == 5
		%Jin, Chiang e Wang
		resp = MAX_mu_X * exp(-0.06035*P -0.0055*S) * (S/(KS + S));
	end
	
	erro = norm(resp - mu_X);
end

function [erro] = OtimizadorODE45(miks_kin, T, X, S, P, t, Pmax, YSX, YPX, esc)
	
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

	options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-4]);
	[~, Y] = ode45(@odefcn, [T(1):2:T(end)], [X0 S0 P0], options);
	
	X1 = Y(:,1);
	S1 = Y(:,2);
	P1 = Y(:,3);
	
	erro = norm(X - X1) + norm(S - S1) + norm(P - P1);
	disp(erro);
end

function [erro] = OtimizadorRK4th(miks_kin, T, X, S, P, t, Pmax, YSX, YPX, esc)
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

function [rsd] = RSD(y, yi, n)
    sum = 0;
	med = 0;
    for i = 1:n
        sum = sum + (y(i) - yi(i))^2;
		med = med + y(i);
    end
	med = med/n;
    aux = sum/n;
	
	rsd = (sqrt(aux)/med)*100;
end

function [Z] = Randomic(y, yi)
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