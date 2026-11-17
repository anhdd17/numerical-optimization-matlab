function [x, fval, status] = solve_lp(c, A, b, Aeq, beq, lb)
% SOLVE_LP  min c'*x  s.t.  A*x <= b,  Aeq*x = beq,  x >= lb
%           bang glpk (GNU Octave) hoac linprog (MATLAB). Bo trong [] neu khong co.
%   status - 1: toi uu; khac 1: vo nghiem / khong bi chan / loi.
  c = c(:);
  n = numel(c);
  if nargin < 6 || isempty(lb), lb = -inf(n, 1); end
  A = reshape(A, [], n);  b = b(:);
  Aeq = reshape(Aeq, [], n);  beq = beq(:);
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    ctype = [repmat('U', 1, size(A, 1)), repmat('S', 1, size(Aeq, 1))];
    [x, fval, err, extra] = glpk(c, [A; Aeq], [b; beq], lb, inf(n, 1), ctype, ...
      repmat('C', 1, n), 1);
    status = double(err == 0 && extra.status == 5);   % 5 = GLP_OPT
  else
    o = optimoptions('linprog', 'Display', 'none');
    [x, fval, flag] = linprog(c, A, b, Aeq, beq, lb, [], o);
    status = double(flag == 1);
  end
end
