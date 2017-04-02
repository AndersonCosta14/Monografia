%% Função com os dados experimentais
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
            disp('Concentrção Inicial inválida');
    end
end