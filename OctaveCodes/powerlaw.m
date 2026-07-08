clear
clc
%% Dados Matheus
#D_km = [5.02, 1.61, 5.0, 2.27, 8.52];  % km
#T_min = [29.65, 7.15, 27.6, 11.9967, 60];  % minutos

%% Dados Rodrigo
D_km = [0.4, 0.8047, 1.0, 1.6093, 3.2187, 5.0, 10.0, 15.0, 16.0934, 20.0, 21.097];
T_min = [1+41/60, 3+37/60, 4+30/60, 7+16/60, 14+56/60, 23+48/60, 52+57/60, 86+42/60, 93+30/60, 117+50/60, 124+57/60];

D = D_km * 1000;
T = T_min * 60;
V = D ./ T
y = log(V);
x = log(T);


A = [ones(length(x), 1), x'];
coef = A \ y';
a = coef(1);
b = coef(2);
S = exp(a)
E = b + 1
F = 1/E

V_ajust = S .* T.^(E-1);
residuos = V - V_ajust;
SS_res = sum(residuos.^2);
SS_tot = sum((V - mean(V)).^2);
R2 = 1 - SS_res/SS_tot

D_5k   = 5000;
D_10k  = 10000;
D_meia = 21097;
D_mara = 42195;

for D_alvo = [D_5k, D_10k, D_meia, D_mara]
    T_alvo = (D_alvo/S)^(1/E);
    fprintf('%.0f m -> %.1f s = %dh%02dmin%02ds\n', ...
        D_alvo, T_alvo, floor(T_alvo/3600), ...
        mod(floor(T_alvo/60), 60), mod(floor(T_alvo), 60))
end
%% Predicoes
T_pred = linspace(60, 8400, 1000);
V_pred = S .* T_pred.^(E-1);
D_pred = S .* T_pred.^E;


%% Plots
figure(1)
% (a) V vs T
subplot(2,2,1)
plot(T_pred, V_pred, 'b-', 'LineWidth', 2)
hold on
scatter(T, V, 30, 'r', 'filled')
xlabel('Tempo de exaustao (s)')
ylabel('Velocidade (m/s)')
title('(a) V vs T - Power Law')
grid on

% (b) log(V) vs log(T)
subplot(2,2,2)
plot(log(T_pred), log(V_pred), 'b-', 'LineWidth', 2)
hold on
scatter(log(T), log(V), 30, 'r', 'filled')
xlabel('log(T)')
ylabel('log(V)')
title('(b) Linearizacao log-log')
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

% (d) V vs 1/T
subplot(2,2,4)
plot(1./T_pred, V_pred, 'b-', 'LineWidth', 2)
hold on
scatter(1./T, V, 30, 'r', 'filled')
xlabel('1/T (1/s)')
ylabel('Velocidade (m/s)')
title('(d) V vs 1/T')
grid on
