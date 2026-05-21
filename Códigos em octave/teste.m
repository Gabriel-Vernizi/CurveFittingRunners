clear
clc
close all

%% Dados Matheus
% D_km = [5.02, 1.61, 5.0, 2.27, 8.52];  % km
% T_min = [29.65, 7.15, 27.6, 11.9967, 60];  % minutos

%% Dados Rodrigo
D_km = [0.4, 0.8047, 1.0, 1.6093, 3.2187, 5.0, 10.0, 15.0, 16.0934, 20.0, 21.097];
T_min = [1+41/60, 3+37/60, 4+30/60, 7+16/60, 14+56/60, 23+48/60, 52+57/60, 86+42/60, 93+30/60, 117+50/60, 124+57/60];

%% Conversao de unidades
D = D_km * 1000;   % metros
T = T_min * 60;    % segundos
V = D ./ T;        % m/s

%% Modelo hipergaussiano
% V(T) = Vinf + A*exp(-(T/tau)^p)
%
% params(1) = Vinf
% params(2) = A
% params(3) = tau
% params(4) = p

hiper_model = @(params, T) params(1) + params(2) .* exp(-(T ./ params(3)).^params(4));

%% Chutes iniciais
Vinf0 = min(V) * 0.85;
A0 = max(V) - Vinf0;
tau0 = median(T);
p0 = 1;

params0 = [Vinf0, A0, tau0, p0];

%% Funcao erro com penalizacao
% A penalizacao evita parametros sem sentido:
% Vinf > 0, A > 0, tau > 0, p > 0

erro = @(params) sum((V - hiper_model(params, T)).^2) + penalizacao(params);

%% Ajuste nao linear
options = optimset('MaxIter', 10000, 'MaxFunEvals', 10000);

params_hiper = fminsearch(erro, params0, options);

Vinf = params_hiper(1);
A_hiper = params_hiper(2);
tau = params_hiper(3);
p = params_hiper(4);

%% Resultados dos parametros
fprintf('\n===== Parametros hipergaussianos =====\n')
fprintf('Vinf = %.4f m/s\n', Vinf)
fprintf('A    = %.4f m/s\n', A_hiper)
fprintf('tau  = %.4f s\n', tau)
fprintf('p    = %.4f\n', p)

%% Valores ajustados nos pontos reais
V_ajust = hiper_model(params_hiper, T);

residuos = V - V_ajust;
SS_res = sum(residuos.^2);
SS_tot = sum((V - mean(V)).^2);
R2 = 1 - SS_res/SS_tot;

RMSE = sqrt(mean(residuos.^2));
MAE = mean(abs(residuos));

fprintf('\n===== Qualidade do ajuste =====\n')
fprintf('R2   = %.6f\n', R2)
fprintf('RMSE = %.6f m/s\n', RMSE)
fprintf('MAE  = %.6f m/s\n', MAE)

%% Predicoes para distancias alvo
D_5k   = 5000;
D_10k  = 10000;
D_meia = 21097;
D_mara = 42195;

fprintf('\n===== Predicoes hipergaussianas =====\n')

for D_alvo = [D_5k, D_10k, D_meia, D_mara]

    func_dist = @(Tteste) Tteste .* hiper_model(params_hiper, Tteste) - D_alvo;

    % intervalo de busca para o tempo
    T_min_busca = 30;
    T_max_busca = 30000;

    T_alvo = fzero(func_dist, [T_min_busca, T_max_busca]);

    fprintf('%.0f m -> %.1f s = %dh%02dmin%02ds\n', ...
        D_alvo, T_alvo, floor(T_alvo/3600), ...
        mod(floor(T_alvo/60), 60), mod(floor(T_alvo), 60))
end

%% Curva para plot
T_pred = linspace(60, 30000, 1500);
V_pred = hiper_model(params_hiper, T_pred);
D_pred = T_pred .* V_pred;

%% Plots
figure(1)

% (a) V vs T
subplot(2,2,1)
plot(T_pred, V_pred, 'b-', 'LineWidth', 2)
hold on
scatter(T, V, 30, 'r', 'filled')
xlabel('Tempo de exaustao (s)')
ylabel('Velocidade (m/s)')
title('(a) V vs T - Hipergaussiana')
grid on

% (b) log(V) vs log(T)
subplot(2,2,2)
plot(log(T_pred), log(V_pred), 'b-', 'LineWidth', 2)
hold on
scatter(log(T), log(V), 30, 'r', 'filled')
xlabel('log(T)')
ylabel('log(V)')
title('(b) Escala log-log')
grid on

% (c) D vs T
subplot(2,2,3)
plot(T_pred, D_pred, 'b-', 'LineWidth', 2)
hold on
scatter(T, D, 30, 'r', 'filled')
xlabel('Tempo de exaustao (s)')
ylabel('Distancia (m)')
title('(c) D vs T')
grid on

% (d) Residuos
subplot(2,2,4)
scatter(T, residuos, 35, 'r', 'filled')
hold on
plot(T_pred, zeros(size(T_pred)), 'k--', 'LineWidth', 1)
xlabel('Tempo de exaustao (s)')
ylabel('Residuo V - V_{ajust}')
title('(d) Residuos')
grid on

%% Tabela comparando dados reais e ajustados
fprintf('\n===== Comparacao ponto a ponto =====\n')
fprintf('Dist(m)\tTempo(s)\tV real\t\tV ajust\t\tResiduo\n')

for i = 1:length(T)
    fprintf('%.0f\t%.1f\t\t%.4f\t\t%.4f\t\t%.4f\n', ...
        D(i), T(i), V(i), V_ajust(i), residuos(i))
end

%% Funcao de penalizacao

