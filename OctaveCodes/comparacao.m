pkg load optim;

% Dados Rodrigo
D_km = [0.4, 0.8047, 1.0, 1.6093, 3.2187, 5.0, 10.0, 15.0, 16.0934, 20.0, 21.097];
T_min = [1+41/60, 3+37/60, 4+30/60, 7+16/60, 14+56/60, 23+48/60, 52+57/60, 86+42/60, 93+30/60, 117+50/60, 124+57/60];

T_seg = T_min * 60;

% Calcular velocidades (km/h)
v_kmh = D_km ./ T_min;

fprintf('\n===============================================\n');
fprintf('DADOS RODRIGO\n');
fprintf('===============================================\n');
fprintf('%-15s %-15s %-15s\n', 'Distância (km)', 'Tempo (min)', 'Vel (km/h)');
fprintf('%-15s %-15s %-15s\n', '---------------', '---------------', '---------------');
for i = 1:length(D_km)
    fprintf('%-15.4f %-15.2f %-15.2f\n', D_km(i), T_min(i), v_kmh(i));
end

% ===================================================================
% MODELO 1: LEI DE POTÊNCIA (RIEGEL)
% ===================================================================
fprintf('\n===============================================\n');
fprintf('RIEGEL (LEI DE POTÊNCIA)\n');
fprintf('===============================================\n');

% Função: v(t) = S * t^(e-1)
riegel = @(T, p) p(1) * T.^(p(2) - 1);

% Palpite inicial
p_riegel_init = [15, 0.08];

% Ajuste
[v_pred_riegel, p_riegel] = leasqr(T_seg, v_kmh, p_riegel_init, riegel);

S_riegel = p_riegel(1);
e_riegel = p_riegel(2);

% Calcular R² e RMSE
residuos_riegel = v_kmh - v_pred_riegel;
ss_res_riegel = sum(residuos_riegel.^2);
ss_tot = sum((v_kmh - mean(v_kmh)).^2);
r2_riegel = 1 - (ss_res_riegel / ss_tot);
rmse_riegel = sqrt(ss_res_riegel / length(v_kmh));

fprintf('S = %.4f km/h\n', S_riegel);
fprintf('e = %.6f\n', e_riegel);
fprintf('R² = %.6f\n', r2_riegel);
fprintf('RMSE = %.6f km/h\n', rmse_riegel);
fprintf('\nInterpretação:\n');
fprintf('  - e = %.4f (e < 0.08 = muito resistente; e > 0.12 = explosivo)\n', e_riegel);
fprintf('  - S = %.4f é coef. escala (v em t=1s)\n', S_riegel);

% ===================================================================
% MODELO 2: SUA HIPERGAUSSIANA
% ===================================================================
fprintf('\n===============================================\n');
fprintf('SEU MODELO (HIPERGAUSSIANA)\n');
fprintf('===============================================\n');

% Função: U(t) = S * g^e / (t^e + g^e)
hypergaussian = @(T, p) p(1) * (p(3).^p(2)) ./ (T.^p(2) + p(3).^p(2));

% Palpite inicial: [S, e, g]
p_hyper_init = [15, 0.08, 100];

% Ajuste
[v_pred_hyper, p_hyper] = leasqr(T_seg, v_kmh, p_hyper_init, hypergaussian);

S_hyper = p_hyper(1);
e_hyper = p_hyper(2);
g_hyper = p_hyper(3);

% Calcular R² e RMSE
residuos_hyper = v_kmh - v_pred_hyper;
ss_res_hyper = sum(residuos_hyper.^2);
r2_hyper = 1 - (ss_res_hyper / ss_tot);
rmse_hyper = sqrt(ss_res_hyper / length(v_kmh));

fprintf('S = %.4f km/h (velocidade inicial)\n', S_hyper);
fprintf('e = %.6f (parâmetro de endurance)\n', e_hyper);
fprintf('g = %.4f segundos (escala temporal)\n', g_hyper);
fprintf('R² = %.6f\n', r2_hyper);
fprintf('RMSE = %.6f km/h\n', rmse_hyper);
fprintf('\nInterpretação:\n');
fprintf('  - S = %.4f km/h é a velocidade inicial (t→0)\n', S_hyper);
fprintf('  - g = %.2f s ≈ %.2f min: tempo de meia-performance\n', g_hyper, g_hyper/60);
fprintf('  - Em t = g: U(g) = S/2 = %.4f km/h\n', S_hyper/2);
fprintf('  - e = %.4f (e ↓ = mais resistente; e ↑ = menos resistente)\n', e_hyper);

% ===================================================================
% COMPARAÇÃO
% ===================================================================
fprintf('\n===============================================\n');
fprintf('COMPARAÇÃO DOS MODELOS\n');
fprintf('===============================================\n');

delta_r2 = r2_hyper - r2_riegel;
delta_rmse = rmse_hyper - rmse_riegel;

fprintf('ΔR² (seu modelo - Riegel) = %+.6f\n', delta_r2);
if delta_r2 > 0.01
    fprintf('  → ✓ Seu modelo melhor\n');
elseif delta_r2 < -0.01
    fprintf('  → ✗ Riegel melhor\n');
else
    fprintf('  → ≈ Equivalentes\n');
end

fprintf('\nΔRMSE = %+.6f km/h\n', delta_rmse);
if delta_rmse < -0.01
    fprintf('  → ✓ Seu modelo melhor\n');
elseif delta_rmse > 0.01
    fprintf('  → ✗ Riegel melhor\n');
else
    fprintf('  → ≈ Equivalentes\n');
end

% AIC (lower is better)
aic_riegel = length(v_kmh) * log(rmse_riegel^2) + 2 * 2;
aic_hyper = length(v_kmh) * log(rmse_hyper^2) + 2 * 3;
fprintf('\nAIC (lower is better):\n');
fprintf('  Riegel: %.2f\n', aic_riegel);
fprintf('  Seu modelo: %.2f\n', aic_hyper);
if aic_hyper < aic_riegel
    fprintf('  → ✓ Seu modelo vence (apesar de 3 parâmetros)\n');
else
    fprintf('  → ✗ Riegel mais eficiente\n');
end

% ===================================================================
% TABELA DE PREDIÇÕES
% ===================================================================
fprintf('\n===============================================\n');
fprintf('PREDIÇÕES vs REAL\n');
fprintf('===============================================\n');
fprintf('%-12s %-12s %-15s %-15s %-12s %-12s\n', 'Tempo', 'Real', 'Riegel', 'Seu modelo', 'Err Rig', 'Err Hyp');
fprintf('%-12s %-12s %-15s %-15s %-12s %-12s\n', '(min)', '(km/h)', '(km/h)', '(km/h)', '(km/h)', '(km/h)');
fprintf('%-12s %-12s %-15s %-15s %-12s %-12s\n', '---', '---', '---', '---', '---', '---');

for i = 1:length(T_seg)
    real = v_kmh(i);
    pred_rig = riegel(T_seg(i), p_riegel);
    pred_hyp = hypergaussian(T_seg(i), p_hyper);
    err_rig = real - pred_rig;
    err_hyp = real - pred_hyp;

    fprintf('%8.2f     %8.3f     %10.3f     %10.3f     %+8.3f     %+8.3f\n', T_min(i), real, pred_rig, pred_hyp, err_rig, err_hyp);
end

% ===================================================================
% PLOTAGEM
% ===================================================================
figure(1);
clf;

% Subplot 1: Velocidade vs Tempo
subplot(1, 2, 1);
plot(T_seg, v_kmh, 'ko', 'markersize', 8, 'linewidth', 2); hold on;

% Plotting range
t_plot = linspace(50, max(T_seg), 200);

% Riegel
v_riegel_plot = riegel(t_plot, p_riegel);
plot(t_plot, v_riegel_plot, 'r-', 'linewidth', 2);

% Hypergaussian
v_hyper_plot = hypergaussian(t_plot, p_hyper);
plot(t_plot, v_hyper_plot, 'b-', 'linewidth', 2);

xlabel('Tempo (segundos)', 'fontsize', 12);
ylabel('Velocidade (km/h)', 'fontsize', 12);
title('Velocidade vs Tempo - Rodrigo', 'fontsize', 13, 'fontweight', 'bold');
legend(sprintf('Dados reais'), sprintf('Riegel (R²=%.4f)', r2_riegel), sprintf('Seu modelo (R²=%.4f)', r2_hyper), 'fontsize', 10);
grid on;
hold off;

% Subplot 2: Resíduos
subplot(1, 2, 2);
plot(T_seg, residuos_riegel, 'rs', 'markersize', 8, 'linewidth', 2); hold on;
plot(T_seg, residuos_hyper, 'bs', 'markersize', 8, 'linewidth', 2);
axhline(y=0, 'color', 'k', 'linestyle', '--', 'linewidth', 1);
xlabel('Tempo (segundos)', 'fontsize', 12);
ylabel('Resíduo (km/h)', 'fontsize', 12);
title('Resíduos (erro = real - predito)', 'fontsize', 13, 'fontweight', 'bold');
legend('Riegel', 'Seu modelo', 'fontsize', 10);
grid on;
hold off;
