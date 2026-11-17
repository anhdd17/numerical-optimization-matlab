function res = lp_dual(c, A, b, display)
% LP_DUAL  Bai toan doi ngau cua quy hoach tuyen tinh dang chinh tac
%
%   (P)  min c'x   s.t.  A x = b,  x >= 0
%   (D)  max b'y   s.t.  A'y <= c          (y tu do)
%
%   Giai ca hai bai toan, kiem tra doi ngau manh (c'x* = b'y*) va dieu kien
%   bu: x_j * (c_j - A_j'y) = 0 voi moi j.
%
%   res = lp_dual(c, A, b)
%   res = lp_dual(c, A, b, display)
%
% Dau ra (struct): .x, .p (nghiem, gia tri cua P); .y, .d (cua D);
%   .slack (= c - A'y, bien bu doi ngau); .compl (= x .* slack); .gap (= p - d)
  if nargin < 4, display = true; end
  c = c(:);  b = b(:);
  [m, n] = size(A);
  [res.x, res.p, sp] = solve_lp(c, [], [], A, b, zeros(n, 1));
  [res.y, negd, sd]  = solve_lp(-b, A', c, [], [], -inf(m, 1));
  res.d     = -negd;
  res.slack = c - A' * res.y;
  res.compl = res.x .* res.slack;
  res.gap   = res.p - res.d;
  res.status = [sp, sd];
  if display
    cl = @(v) v .* (abs(v) > 1e-9) + 0;        % bo nhieu lam tron khi in
    fprintf('  (P): x* = %s,  p* = %g\n', mat2str(cl(res.x)', 4), res.p);
    fprintf('  (D): y* = %s,  d* = %g\n', mat2str(cl(res.y)', 4), res.d);
    fprintf('  Khe ho p* - d* = %.2e\n', res.gap);
    fprintf('  Bien bu doi ngau c - A''y* = %s\n', mat2str(cl(res.slack)', 4));
    fprintf('  Dieu kien bu x_j*(c - A''y)_j = %s\n', mat2str(cl(res.compl)', 4));
  end
end
