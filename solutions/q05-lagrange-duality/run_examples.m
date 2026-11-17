% RUN_EXAMPLES  Cac vi du minh hoa ham doi ngau Lagrange va cac dinh ly doi ngau (cau 5).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;

%% Vi du 1: quy hoach toan phuong (bai toan loi, Slater thoa) - doi ngau manh
fprintf('==================== VI DU 1 ====================\n');
fprintf('(a) min x1^2 + x2^2  s.t.  x1 + x2 >= 1   (viet lai: -x1 - x2 <= -1)\n');
Q = 2 * eye(2);  c = [0; 0];
r1 = dual_qp(Q, c, [-1 -1], -1);
fprintf('  q(lambda) = lambda - lambda^2/2;  kiem tra q(0.5) = %g, q(1) = %g, q(2) = %g\n', ...
  r1.q(0.5, []), r1.q(1, []), r1.q(2, []));
fprintf('  (P): x* = %s, p* = %g\n', mat2str(r1.xP', 4), r1.p);
fprintf('  (D): lambda* = %g, d* = %g, x(lambda*) = %s\n', r1.lambda, r1.d, mat2str(r1.xD', 4));
fprintf('  Khe ho p* - d* = %.2e\n\n', r1.gap);

fprintf('(b) min x1^2 + x2^2  s.t.  x1 + 2x2 = 5   (rang buoc dang thuc, mu tu do)\n');
r1b = dual_qp(Q, c, [], [], [1 2], 5);
fprintf('  q(mu) = -5mu^2/4 - 5mu;  kiem tra q(-2) = %g\n', r1b.q([], -2));
fprintf('  (P): x* = %s, p* = %g\n', mat2str(r1b.xP', 4), r1b.p);
fprintf('  (D): mu* = %g, d* = %g, x(mu*) = %s\n', r1b.mu, r1b.d, mat2str(r1b.xD', 4));
fprintf('  Khe ho p* - d* = %.2e\n\n', r1b.gap);

%% Vi du 2: quy hoach tuyen tinh (vi du 1 cua cau 2 va cau 3)
%   min 4x1 + x2  s.t.  3x1 + x2 - x3 = 3,  4x1 + 3x2 - x4 = 6,  x1 + 2x2 + x5 = 4,  x >= 0
fprintf('==================== VI DU 2 ====================\n');
c2 = [4; 1; 0; 0; 0];
A2 = [3 1 -1  0 0;
      4 3  0 -1 0;
      1 2  0  0 1];
b2 = [3; 6; 4];
r2 = lp_dual(c2, A2, b2);
fprintf('  Kiem tra rang buoc doi ngau: A''y* <= c  -> %s\n\n', ...
  mat2str(all(A2' * r2.y <= c2 + 1e-9)));

%% Vi du 3: bai toan khong loi - khe ho doi ngau duong
%   min -x^2  s.t.  x - 1 <= 0,  x in X = [0, 2]
%   q(lambda) = min_{x in [0,2]} -x^2 + lambda(x - 1) = min(-lambda, lambda - 4)
fprintf('==================== VI DU 3 ====================\n');
xg = linspace(0, 2, 2001);                    % tinh q bang luoi tren X
L3 = linspace(0, 4, 4001);
q3 = zeros(size(L3));
for k = 1:numel(L3)
  q3(k) = min(-xg.^2 + L3(k) * (xg - 1));
end
[d3, k3] = max(q3);
p3 = -1;                                      % x* = 1
fprintf('  p* = %g (x* = 1)\n', p3);
fprintf('  d* = %g tai lambda* = %g (so: tinh tren luoi; giai tich: d* = -2, lambda* = 2)\n', ...
  d3, L3(k3));
fprintf('  Khe ho doi ngau p* - d* = %g > 0\n', p3 - d3);
fprintf('  Sai lech giua q tren luoi va cong thuc min(-l, l-4): %.2e\n\n', ...
  max(abs(q3 - min(-L3, L3 - 4))));

%% Vi du 4: bai toan loi khong thoa Slater (vi du 6 cua cau 4)
%   min x  s.t.  x^2 <= 0.   q(lambda) = -1/(4 lambda) voi lambda > 0, q(0) = -inf
fprintf('==================== VI DU 4 ====================\n');
for lam = [1 10 100 1e4]
  % so: x(lambda) = -1/(2 lambda) la diem cuc tieu cua x + lambda x^2
  xl = -1 / (2 * lam);
  fprintf('  lambda = %-6g: q = %-12g (= -1/(4 lambda) = %g)\n', lam, xl + lam * xl^2, -1 / (4 * lam));
end
fprintf('  p* = 0, sup q = 0 nhung khong dat duoc: bai toan doi ngau khong co nghiem.\n\n');

%% Hinh: ham doi ngau cua vi du 1(a), 3, 4
if ~exist('figures', 'dir'), mkdir('figures'); end
fig = figure('visible', 'off', 'Position', [100 100 1500 460]);

subplot(1, 3, 1); hold on; grid on; box on;
l = linspace(0, 2.5, 200);
h1 = plot(l, l - l.^2 / 2, 'b-', 'LineWidth', 1.8);
h2 = plot([0 2.5], [0.5 0.5], 'r--', 'LineWidth', 1.2);
plot(1, 0.5, 'ko', 'MarkerFaceColor', 'k');
text(0.15, 0.66, '\lambda^* = 1,  d^* = p^* = 1/2', 'FontSize', 10);
axis([0 2.5 -0.6 0.8]);
xlabel('\lambda'); ylabel('q(\lambda)');
legend([h1, h2], {'q(\lambda) = \lambda - \lambda^2/2', 'p^* = 1/2'}, 'Location', 'southwest');
title('Vi du 1(a): khong co khe ho');

subplot(1, 3, 2); hold on; grid on; box on;
h1 = plot(L3, q3, 'b-', 'LineWidth', 1.8);
h2 = plot([0 4], [p3 p3], 'r--', 'LineWidth', 1.2);
plot(2, -2, 'ko', 'MarkerFaceColor', 'k');
plot([2 2], [-2 -1], 'k:', 'LineWidth', 1.5);
text(2.12, -1.5, 'khe ho = 1', 'FontSize', 10);
text(1.15, -1.8, 'd^* = -2', 'FontSize', 10);
axis([0 4 -4.2 0.3]);
xlabel('\lambda'); ylabel('q(\lambda)');
legend([h1, h2], {'q(\lambda) = min(-\lambda, \lambda - 4)', 'p^* = -1'}, 'Location', 'north');
title('Vi du 3: khe ho doi ngau duong');

subplot(1, 3, 3); hold on; grid on; box on;
l = linspace(0.05, 5, 300);
h1 = plot(l, -1 ./ (4 * l), 'b-', 'LineWidth', 1.8);
h2 = plot([0 5], [0 0], 'r--', 'LineWidth', 1.2);
axis([0 5 -2.6 0.3]);
xlabel('\lambda'); ylabel('q(\lambda)');
legend([h1, h2], {'q(\lambda) = -1/(4\lambda)', 'p^* = 0 = sup q (khong dat)'}, ...
  'Location', 'southeast');
title('Vi du 4: d^* = p^* nhung khong dat');

print(fig, '-dpng', '-r130', fullfile('figures', 'dual_functions.png'));
close(fig);
fprintf('Da luu hinh figures/dual_functions.png\n');
