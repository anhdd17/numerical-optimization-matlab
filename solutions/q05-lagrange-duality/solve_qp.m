function [x, fval, lambda] = solve_qp(H, f, A, b, Aeq, beq, lb)
% SOLVE_QP  min 0.5*x'*H*x + f'*x  s.t.  A*x <= b,  Aeq*x = beq,  x >= lb
%           bang qp (GNU Octave) hoac quadprog (MATLAB). Bo trong [] neu khong co.
%   lambda - (Octave) nhan tu Lagrange cua cac rang buoc, theo thu tu cua qp.
  n = numel(f);
  if nargin < 7 || isempty(lb), lb = -inf(n, 1); end
  if nargin < 6, Aeq = []; beq = []; end
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    Ain = A; bin = b;
    if isempty(Ain), Ain = []; bin = []; end
    [x, fval, ~, lambda] = qp(zeros(n, 1), H, f(:), Aeq, beq, lb, [], [], Ain, bin);
  else
    o = optimoptions('quadprog', 'Display', 'none');
    [x, fval, ~, ~, l] = quadprog(H, f(:), A, b, Aeq, beq, lb, [], [], o);
    lambda = [l.eqlin; l.lower; l.ineqlin];
  end
end
