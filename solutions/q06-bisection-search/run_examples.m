% RUN_EXAMPLES  Cac vi du minh hoa giai thuat chia doi tim cuc tieu ham mot bien (cau 6).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac 4 unimodal tren [0, 2]
%   f(x) = x^4 - 14x^3 + 60x^2 - 70x
fprintf('==================== VI DU 1 ====================\n');
f1  = @(x) x.^4 - 14*x.^3 + 60*x.^2 - 70*x;
df1 = @(x) 4*x.^3 - 42*x.^2 + 120*x - 70;
fprintf('--- Chia doi dung gia tri ham (tol = 1e-5, delta = 2.5e-6) ---\n');
[x1, fx1, o1] = bisection_min(f1, 0, 2, struct('tol', 1e-5, 'display', true));
fprintf('x* = %.8f,  f* = %.8f,  so buoc = %d,  so lan tinh f = %d\n\n', ...
  x1, fx1, o1.iterations, o1.nfev);
[xr, fr] = fminbnd(f1, 0, 2, optimset('TolX', 1e-10));
fprintf('Doi chieu fminbnd: x* = %.8f, f* = %.8f, |x - x_ref| = %.1e\n\n', xr, fr, abs(x1 - xr));

%% Vi du 2: so sanh hai bien the tren f(x) = e^{-x} + x^2, [0, 1]
fprintf('==================== VI DU 2 ====================\n');
f2  = @(x) exp(-x) + x.^2;
df2 = @(x) -exp(-x) + 2*x;
tol = 1e-6;
[xa, ~, oa] = bisection_min(f2, 0, 1, struct('tol', tol));
[xb, ~, ob] = bisection_min(f2, 0, 1, struct('tol', tol, 'method', 'derivative', 'df', df2));
xr2 = fzero(df2, [0 1]);                      % nghiem cua f'(x) = 0
delta = tol / 4;
kA = ceil(log2((1 - 2*delta) / (tol - 2*delta)));
kB = ceil(log2(1 / tol));
fprintf('%-22s %12s %8s %8s %10s %12s\n', 'Bien the', 'x*', 'buoc', 'ly thuyet', 'tinh f/f''', '|x - x_ref|');
fprintf('%-22s %12.8f %8d %8d %10d %12.1e\n', 'gia tri ham', xa, oa.iterations, kA, oa.nfev, abs(xa - xr2));
fprintf('%-22s %12.8f %8d %8d %10d %12.1e\n', 'dao ham', xb, ob.iterations, kB, ob.ndfev, abs(xb - xr2));
fprintf('Nghiem tham chieu (fzero cua f''): x* = %.8f\n\n', xr2);

%% Vi du 3: ham khong unimodal - chia doi chi tim duoc cuc tieu dia phuong
%   f(x) = sin(x) + sin(10x/3) tren [2.7, 7.5]
fprintf('==================== VI DU 3 ====================\n');
f3 = @(x) sin(x) + sin(10*x/3);
[x3, fx3] = bisection_min(f3, 2.7, 7.5);
fprintf('Chia doi tren ca doan [2.7, 7.5]: x = %.6f, f = %.6f (cuc tieu dia phuong)\n', x3, fx3);
% Khac phuc: quet luoi tho de khoanh vung, roi chia doi tren khoang quanh diem luoi tot nhat
N  = 25;
xg = linspace(2.7, 7.5, N);
[~, i] = min(f3(xg));
lo = xg(max(i - 1, 1));  hi = xg(min(i + 1, N));
[x3g, fx3g] = bisection_min(f3, lo, hi);
fprintf('Quet luoi %d diem -> khoang [%.4f, %.4f] -> chia doi: x = %.6f, f = %.6f\n', ...
  N, lo, hi, x3g, fx3g);
xs = linspace(2.7, 7.5, 200001);
[fm, im] = min(f3(xs));
fprintf('Tham chieu (luoi 200001 diem): x = %.6f, f = %.6f\n\n', xs(im), fm);

%% Hinh 1: vi du 1 - ham va cac khoang qua 6 buoc dau
fig = figure('visible', 'off', 'Position', [100 100 900 560]);
subplot(2, 1, 1); hold on; grid on; box on;
xx = linspace(0, 2, 400);
plot(xx, f1(xx), 'b-', 'LineWidth', 1.6);
plot(x1, fx1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(x1 + 0.07, fx1 + 7, sprintf('x^* = %.4f', x1), 'FontSize', 10);
plot([x1 x1 + 0.06], [fx1 + 0.8, fx1 + 6], 'k-');
xlim([0 2]);
ylabel('f(x)');
title('Vi du 1: f(x) = x^4 - 14x^3 + 60x^2 - 70x va cac khoang [a_k, b_k]');
subplot(2, 1, 2); hold on; grid on; box on;
K = 6;
for k = 0:K
  plot([o1.a(k+1) o1.b(k+1)], [k k], 'k-', 'LineWidth', 4);
  if k == 0
    text(1.55, 0.45, sprintf('L_%d = %.4f', k, o1.b(k+1) - o1.a(k+1)), 'FontSize', 9);
  else
    text(o1.b(k+1) + 0.03, k, sprintf('L_%d = %.4f', k, o1.b(k+1) - o1.a(k+1)), 'FontSize', 9);
  end
end
plot([x1 x1], [-0.5 K + 0.5], 'r--');
set(gca, 'YDir', 'reverse');
axis([0 2 -0.5 K + 0.5]);
xlabel('x'); ylabel('buoc k');
print(fig, '-dpng', '-r150', fullfile('figures', 'example1_intervals.png'));
close(fig);

%% Hinh 2: toc do hoi tu - do dai khoang theo so lan tinh ham (vi du 2)
fig = figure('visible', 'off');
LA = oa.b - oa.a;   LB = ob.b - ob.a;
semilogy(0:2:2*(numel(LA) - 1), LA, 'bo-', 'MarkerSize', 4); hold on; grid on;
semilogy(0:numel(LB) - 1, LB, 'rs-', 'MarkerSize', 4);
kk = 0:numel(LA) - 1;
semilogy(2*kk, (1 - 2*delta) ./ 2.^kk + 2*delta, 'k--');
semilogy([0 45], [tol tol], 'k:');
text(30, 2.2e-6, 'tol = 10^{-6}', 'FontSize', 10);
xlim([0 45]);
xlabel('so lan tinh f (hoac f'')'); ylabel('do dai khoang L');
legend('Gia tri ham (2 lan tinh f / buoc)', 'Dao ham (1 lan tinh f'' / buoc)', ...
  'Ly thuyet: (L_0 - 2\delta)/2^k + 2\delta', 'Location', 'northeast');
title('Vi du 2: do dai khoang giam theo cap so nhan cong boi 1/2');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_convergence.png'));
close(fig);

%% Hinh 3: vi du 3 - ham nhieu cuc tieu
fig = figure('visible', 'off');
hold on; grid on; box on;
xx = linspace(2.7, 7.5, 600);
plot(xx, f3(xx), 'b-', 'LineWidth', 1.6);
plot(xg, f3(xg), 'k.', 'MarkerSize', 12);
plot(x3, fx3, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(x3g, fx3g, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(x3 - 0.25, fx3 - 0.25, 'chia doi tren [2.7, 7.5]', 'FontSize', 9, 'Color', 'r');
text(x3g + 0.15, fx3g - 0.05, 'quet luoi + chia doi', 'FontSize', 9);
axis([2.7 7.5 -2.2 2.2]);
xlabel('x'); ylabel('f(x)');
legend('f(x) = sin x + sin(10x/3)', 'Luoi tho 25 diem', 'Cuc tieu dia phuong', ...
  'Cuc tieu toan cuc', 'Location', 'north');
title('Vi du 3: ham khong unimodal');
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_multimodal.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');
