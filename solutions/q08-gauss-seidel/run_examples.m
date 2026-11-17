% RUN_EXAMPLES  Cac vi du minh hoa giai thuat Gauss-Seidel (toa do tuan hoan) (cau 8).
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> run_examples
clear; clc;
addpath(fullfile('..', 'q07-golden-section-search'));   % golden_section.m
if ~exist('figures', 'dir'), mkdir('figures'); end

%% Vi du 1: ham bac hai - trung voi phep lap Gauss-Seidel giai A x = b
%   f(x) = 1/2 x'Ax - b'x,  A = [2 1; 1 3],  b = [1; 2],  x0 = (0, 0)
fprintf('==================== VI DU 1 ====================\n');
A1 = [2 1; 1 3];  b1 = [1; 2];
f1 = @(x) 0.5 * x' * A1 * x - b1' * x;
% (a) tim kiem mot chieu bang lat cat vang
[xa, fa, oa] = gauss_seidel_min(f1, [0; 0], struct('tol', 1e-8));
% (b) cuc tieu chinh xac theo x_i: df/dx_i = 0  =>  x_i = (b_i - sum_{j~=i} a_ij x_j) / a_ii
cm = @(x, i) (b1(i) - A1(i, :) * x + A1(i, i) * x(i)) / A1(i, i);
[xb, fb, ob] = gauss_seidel_min(f1, [0; 0], struct('tol', 1e-8, 'coordMin', cm));
% (c) phep lap Gauss-Seidel cho he tuyen tinh: x^{k+1} = (D + L)^{-1} (b - U x^k)
DL = tril(A1);  U = triu(A1, 1);
xc = [0; 0];  Xc = xc;
for k = 1:ob.iterations
  xc = DL \ (b1 - U * xc);
  Xc(:, end+1) = xc;
end
fprintf('%3s | %-24s | %-24s | %-24s\n', 'k', '(a) lat cat vang', '(b) cong thuc', '(c) GS cho Ax = b');
for k = 0:4
  fprintf('%3d | %-24s | %-24s | %-24s\n', k, mat2str(oa.X(:, k+1)', 6), ...
    mat2str(ob.X(:, k+1)', 6), mat2str(Xc(:, k+1)', 6));
end
fprintf('Sai lech lon nhat (a)-(c) tren 5 vong dau: %.1e;  (b)-(c): %.1e\n', ...
  max(max(abs(oa.X(:, 1:5) - Xc(:, 1:5)))), max(max(abs(ob.X - Xc))));
rho1 = max(abs(eig(DL \ U)));
fprintf('Nghiem A\\b = %s; (b) dung sau %d vong; ban kinh pho cua ma tran lap = %.4f\n\n', ...
  mat2str((A1 \ b1)', 6), ob.iterations, rho1);

%% Vi du 2: ham khong bac hai (vi du 3 cua cau 7) - thung lung hep
%   F(x) = (x1 - 2)^4 + (x1 - 2x2)^2,  x0 = (0, 3),  x* = (2, 1)
fprintf('==================== VI DU 2 ====================\n');
F2 = @(x) (x(1) - 2)^4 + (x(1) - 2*x(2))^2;
[~, ~, o2s] = gauss_seidel_min(F2, [0; 3], struct('display', true, 'maxIter', 6));
[x2, f2, o2] = gauss_seidel_min(F2, [0; 3], struct('maxIter', 500));
fprintf('Sau %d vong: x = %s, F = %.3e, ||x - x*|| = %.3e, %d lan tinh F\n', ...
  o2.iterations, mat2str(x2', 6), f2, norm(x2 - [2; 1]), o2.nfev);
for K = [10 50 100 500]
  fprintf('  vong %3d: ||x - x*|| = %.3e, F = %.3e\n', K, norm(o2.X(:, K+1) - [2; 1]), o2.F(K+1));
end
fprintf('\n');

%% Vi du 3: n = 10 bien - ham bac hai voi ma tran ba duong cheo
%   f(x) = 1/2 x'Ax - b'x,  A = tridiag(-1, 2, -1),  b = (1, ..., 1)
fprintf('==================== VI DU 3 ====================\n');
n = 10;
A3 = 2 * eye(n) - diag(ones(n - 1, 1), 1) - diag(ones(n - 1, 1), -1);
b3 = ones(n, 1);
f3 = @(x) 0.5 * x' * A3 * x - b3' * x;
xs3 = A3 \ b3;                                % nghiem: x_i = i(n + 1 - i)/2
cm3 = @(x, i) (b3(i) - A3(i, :) * x + A3(i, i) * x(i)) / A3(i, i);
% Chi dung tieu chuan ||dx|| <= tol (ftol = 0): gan cuc tieu f giam theo binh phuong sai so,
% nen tieu chuan theo f se dung qua som khi rho gan 1 (xem them dong in o duoi).
[x3b, ~, o3b] = gauss_seidel_min(f3, zeros(n, 1), struct('tol', 1e-8, 'ftol', 0, 'coordMin', cm3, 'maxIter', 2000));
[x3f, ~, o3f] = gauss_seidel_min(f3, zeros(n, 1), struct('tol', 1e-8, 'coordMin', cm3, 'maxIter', 2000));
[x3a, ~, o3a] = gauss_seidel_min(f3, zeros(n, 1), struct('tol', 1e-8, 'ftol', 0, 'maxIter', 2000));
rho3 = max(abs(eig(tril(A3) \ triu(A3, 1))));
fprintf('Nghiem chinh xac x_i = i(n+1-i)/2 = %s\n', mat2str(round(xs3')));
fprintf('Cong thuc toa do:  %d vong, ||x - x*|| = %.2e\n', o3b.iterations, norm(x3b - xs3));
fprintf('Lat cat vang:      %d vong, ||x - x*|| = %.2e, %d lan tinh f\n', ...
  o3a.iterations, norm(x3a - xs3), o3a.nfev);
fprintf('Neu them tieu chuan dung theo f (ftol = 1e-12): dung sau %d vong, ||x - x*|| = %.2e\n', ...
  o3f.iterations, norm(x3f - xs3));
fprintf('Ban kinh pho rho = %.4f (= cos^2(pi/(n+1)) = %.4f)\n\n', rho3, cos(pi / (n + 1))^2);

%% Vi du 4: ham khong tron - Gauss-Seidel bi ket
%   f(x) = 2|x1 - x2| + (x1 + x2 - 3)^2,  x0 = (0, 0),  x* = (1.5, 1.5), f* = 0
fprintf('==================== VI DU 4 ====================\n');
f4 = @(x) 2 * abs(x(1) - x(2)) + (x(1) + x(2) - 3)^2;
[x4, fx4, o4] = gauss_seidel_min(f4, [0; 0], struct('display', true));
fprintf('Dung tai x = %s, f = %.6f; cuc tieu that x* = (1.5, 1.5), f* = 0\n', mat2str(x4', 6), fx4);
fprintf('Tai (2,2): f(2+t, 2) va f(2, 2+t) voi t = -0.1, 0.1: %s\n\n', ...
  mat2str([f4([1.9; 2]), f4([2.1; 2]), f4([2; 1.9]), f4([2; 2.1])], 4));

%% Hinh 1: vi du 1 va vi du 4 - duong zic-zac tren duong muc
fig = figure('visible', 'off', 'Position', [100 100 1100 460]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.4, 1, 200), linspace(-0.2, 1.2, 200));
Z = 0.5 * (A1(1,1) * X1.^2 + 2 * A1(1,2) * X1 .* X2 + A1(2,2) * X2.^2) - b1(1) * X1 - b1(2) * X2;
contour(X1, X2, Z, f1(A1 \ b1) + [0.01 0.05 0.15 0.3 0.5 0.8], 'LineColor', [0.6 0.6 0.6]);
P = ob.path;
plot(P(1, :), P(2, :), 'b.-', 'LineWidth', 1.4, 'MarkerSize', 12);
plot(P(1, 1), P(2, 1), 'ko', 'MarkerFaceColor', 'k');
plot(xb(1), xb(2), 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(0.03, -0.1, 'x_0', 'FontSize', 11);
text(0.25, 0.68, 'x^* = A^{-1}b', 'FontSize', 10, 'BackgroundColor', 'w');
axis equal; axis([-0.4 1 -0.2 1.2]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 1: f bac hai, cac buoc doc truc toa do');
subplot(1, 2, 2); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.3, 3, 300), linspace(-0.3, 3, 300));
Z = 2 * abs(X1 - X2) + (X1 + X2 - 3).^2;
contour(X1, X2, Z, [0.25 0.5 1 1.5 2.5 4 6 9], 'LineColor', [0.6 0.6 0.6]);
plot([-0.3 3], [-0.3 3], 'k:');
P = o4.path;
plot(P(1, :), P(2, :), 'r.-', 'LineWidth', 1.4, 'MarkerSize', 12);
plot(2, 2, 'rv', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(1.5, 1.5, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(2.1, 1.85, 'ket tai (2, 2), f = 1', 'Color', 'r', 'FontSize', 10, 'BackgroundColor', 'w');
plot([1.5 0.95], [1.5 2.3], 'k-');
text(0.2, 2.42, 'x^* = (1.5, 1.5), f^* = 0', 'FontSize', 10, 'BackgroundColor', 'w');
text(0.08, 0.12, 'x_0', 'FontSize', 11);
axis equal; axis([-0.3 3 -0.3 3]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 4: f khong tron, Gauss-Seidel bi ket');
print(fig, '-dpng', '-r150', fullfile('figures', 'examples1_4_paths.png'));
close(fig);

%% Hinh 2: vi du 2 - duong di va toc do hoi tu
fig = figure('visible', 'off', 'Position', [100 100 1100 460]);
subplot(1, 2, 1); hold on; grid on; box on;
[X1, X2] = meshgrid(linspace(-0.5, 3.5, 300), linspace(0, 3.5, 300));
Z = (X1 - 2).^4 + (X1 - 2 * X2).^2;
contour(X1, X2, Z, [0.01 0.1 0.37 1 3 8 16 30 52], 'LineColor', [0.6 0.6 0.6]);
P = o2.path(:, 1:2 * 20 + 1);
plot(P(1, :), P(2, :), 'b.-', 'LineWidth', 1.2, 'MarkerSize', 10);
plot(0, 3, 'ko', 'MarkerFaceColor', 'k');
plot(2, 1, 'kp', 'MarkerSize', 13, 'MarkerFaceColor', [1 0.8 0]);
text(0.1, 3.2, 'x_0 = (0, 3)', 'FontSize', 10);
text(1.3, 0.55, 'x^* = (2, 1)', 'FontSize', 10, 'BackgroundColor', 'w');
axis equal; axis([-0.5 3.5 0 3.5]);
xlabel('x_1'); ylabel('x_2');
title('Vi du 2: 20 vong dau');
subplot(1, 2, 2);
K = o2.iterations;
err = sqrt(sum((o2.X - [2; 1]).^2, 1));
loglog(1:K, err(2:end), 'b-', 'LineWidth', 1.6); hold on; grid on;
loglog(1:K, o2.F(2:end), 'r-', 'LineWidth', 1.6);
xlabel('vong lap k'); ylabel('gia tri');
legend('||x_k - x^*||', 'F(x_k) - F^*', 'Location', 'southwest');
title('Vi du 2: hoi tu cham (thang log-log)');
print(fig, '-dpng', '-r150', fullfile('figures', 'example2_valley.png'));
close(fig);

%% Hinh 3: vi du 3 - sai so theo vong lap, so voi rho^k
fig = figure('visible', 'off');
e = sqrt(sum((o3b.X - xs3).^2, 1));
semilogy(0:numel(e) - 1, e, 'b-', 'LineWidth', 1.6); hold on; grid on;
kk = 0:numel(e) - 1;
semilogy(kk, e(1) * rho3.^kk, 'k--');
xlabel('vong lap k'); ylabel('||x_k - x^*||');
legend('Gauss-Seidel (n = 10)', sprintf('||x_0 - x^*|| \\rho^k, \\rho = %.3f', rho3), ...
  'Location', 'northeast');
title('Vi du 3: hoi tu tuyen tinh voi he so \rho');
print(fig, '-dpng', '-r150', fullfile('figures', 'example3_convergence.png'));
close(fig);
fprintf('Da luu hinh vao thu muc figures/\n');

