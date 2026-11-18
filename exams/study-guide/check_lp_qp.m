% CHECK_LP_QP  Tinh dap an cac cau hoi mau trong topics/lp-and-qp.md
% (cau 1, 2, 3, 18 cua bai tap lon). Moi so trong tai lieu lay tu output nay.
% Chay duoc tren ca MATLAB va GNU Octave:
%   >> check_lp_qp
clear; clc;
exam_setup;
quiet = struct('display', false);

%% Cau 1 - mau 1.1: max F = 3x1 + 5x2, x1 <= 4, 2x2 <= 12, 3x1 + 2x2 <= 18
fprintf('==================== MAU 1.1 ====================\n');
c11 = [-3; -5; 0; 0; 0];                       % min -F
A11 = [1 0 1 0 0; 0 2 0 1 0; 3 2 0 0 1];  b11 = [4; 12; 18];
[x, fv, ef, o] = simplex(c11, A11, b11, [3 4 5], struct('display', true));
fprintf('Don hinh: x* = %s, F* = %g, %d buoc, exitflag = %d\n', mat2str(x' + 0, 6), -fv, o.iterations, ef);
[xl, fl] = solve_lp(c11, [], [], A11, b11, zeros(5, 1));
fprintf('Solver: x = %s, F = %g\n\n', mat2str(xl', 6), -fl);

%% Cau 1 - mau 1.2: min f = -x1 - 2x2, x1 - x2 + x3 = 2, -x1 + x2 + x4 = 3
fprintf('==================== MAU 1.2 ====================\n');
c12 = [-1; -2; 0; 0];  A12 = [1 -1 1 0; -1 1 0 1];  b12 = [2; 3];
[x, fv, ef, o] = simplex(c12, A12, b12, [3 4], struct('display', true));
fprintf('exitflag = %d, huong d = %s, c''d = %g, A*d = %s\n', ef, mat2str(o.direction' + 0), c12' * o.direction, mat2str((A12 * o.direction)'));
[xl, fl, st] = solve_lp(c12, [], [], A12, b12, zeros(4, 1));
fprintf('Solver status = %d (khac 1: khong co nghiem toi uu)\n\n', st);

%% Cau 1 - mau 1.3: vi du 1 cau 1, phuong an toi uu co suy bien khong
fprintf('==================== MAU 1.3 ====================\n');
c13 = [-3; -2; 0; 0; 0];  A13 = [1 1 1 0 0; 1 3 0 1 0; 1 0 0 0 1];  b13 = [4; 6; 3];
[x, fv, ef, o] = simplex(c13, A13, b13, [], quiet);
fprintf('x* = %s, f* = %g, co so cuoi = %s, x_B = %s\n\n', mat2str(x' + 0, 6), fv, mat2str(o.basis), mat2str(x(o.basis)' + 0, 6));

%% Cau 2 - mau 2.1: min 2x1 + 3x2, x1 + x2 >= 4, x1 + 3x2 >= 6
fprintf('==================== MAU 2.1 ====================\n');
c21 = [2; 3; 0; 0];  A21 = [1 1 -1 0; 1 3 0 -1];  b21 = [4; 6];
[x, fv, ef, o] = two_phase_simplex(c21, A21, b21, struct('display', true));
fprintf('Hai pha: x* = %s, f* = %g, bien gia = %s, w* = %g, buoc pha 1 = %d, pha 2 = %d\n', ...
  mat2str(x' + 0, 6), fv, mat2str(o.artificial), o.phase1.w, o.phase1.iterations, o.phase2.iterations);
[xl, fl] = solve_lp(c21, [], [], A21, b21, zeros(4, 1));
fprintf('Solver: x = %s, f = %g\n', mat2str(xl', 6), fl);
y21 = A21(:, 1:2)' \ c21(1:2);
fprintf('Bien doi ngau tai co so {x1, x2}: y = %s\n\n', mat2str(y21', 6));

%% Cau 2 - mau 2.2: min x1 + x2, x1 + x2 >= 5, x1 + 2x2 <= 4 (vo nghiem)
fprintf('==================== MAU 2.2 ====================\n');
c22 = [1; 1; 0; 0];  A22 = [1 1 -1 0; 1 2 0 1];  b22 = [5; 4];
[x, fv, ef, o] = two_phase_simplex(c22, A22, b22, struct('display', true));
fprintf('Hai pha: exitflag = %d, w* = %g, bien gia = %s\n', ef, o.phase1.w, mat2str(o.artificial));
[xl, fl, st] = solve_lp(c22, [], [], A22, b22, zeros(4, 1));
fprintf('Solver status = %d\n\n', st);

%% Cau 3 - mau 3.1, 3.2: bai toan M cho bai toan cua mau 2.1 / 2.2
fprintf('==================== MAU 3.1 / 3.2 ====================\n');
[x, fv, ef, o] = big_m_simplex(c21, A21, b21, struct('display', true));
fprintf('Bai toan M (M hinh thuc): x* = %s, f* = %g, exitflag = %d, %d buoc\n', mat2str(x' + 0, 6), fv, ef, o.iterations);
for M = [1 1.4 1.6 2 100]
  [x, fv, ef, o] = big_m_simplex(c21, A21, b21, struct('display', false, 'M', M));
  fprintf('  M = %5.1f: exitflag = %2d, x = %s, bien gia = %s, f = %g\n', M, ef, mat2str(x' + 0, 6), mat2str(o.xa' + 0, 6), fv);
end
[x, fv, ef, o] = big_m_simplex(c22, A22, b22, struct('display', false));
fprintf('Mau 2.2 bang bai toan M: exitflag = %d, bien gia = %s\n\n', ef, mat2str(o.xa' + 0, 6));

%% Cau 18 - mau 18.1: min x1^2 + x2^2 - 4x1 - 2x2, x1 + x2 + x3 = 2, x >= 0
fprintf('==================== MAU 18.1 ====================\n');
c18 = [-4; -2; 0];  D18 = diag([2 2 0]);  A18 = [1 1 1];  b18 = 2;
[x, fv, ef, o] = beale_qp(c18, D18, A18, b18, struct('display', true, 'B0', 3));
fprintf('Beale: x* = %s, f* = %g, exitflag = %d, %d buoc, nhan tu df/dx = %s\n', ...
  mat2str(x' + 0, 6), fv, ef, o.iterations, mat2str(o.lambda' + 0, 6));
[xq, fq] = solve_qp(diag([2 2]), [-4; -2], [1 1], 2, [], [], zeros(2, 1));
fprintf('Solver qp (2 bien, x1 + x2 <= 2): x = %s, f = %g\n', mat2str(xq', 6), fq);
% KKT tay: 2x1 - 4 + mu = 0, 2x2 - 2 + mu = 0, x1 + x2 = 2
fprintf('KKT tay: mu = %g\n\n', 4 - 2 * x(1));

%% Cau 18 - mau 18.2: vi du 1 cau 18 (bai toan cua cau 15, 17)
fprintf('==================== MAU 18.2 ====================\n');
c182 = [-4; -6; 0; 0];  D182 = [4 -2 0 0; -2 4 0 0; 0 0 0 0; 0 0 0 0];
A182 = [1 1 1 0; 1 5 0 1];  b182 = [2; 5];
[x, fv, ef, o] = beale_qp(c182, D182, A182, b182, struct('display', true, 'B0', [3 4]));
fprintf('Beale: x* = %s, f* = %.6f (= -222/31 = %.6f), %d buoc, so bien tu do = %d\n', ...
  mat2str(x' + 0, 6), fv, -222/31, o.iterations, o.nFree);
