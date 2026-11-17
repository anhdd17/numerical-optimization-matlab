function res = dual_qp(Q, c, A, b, Aeq, beq)
% DUAL_QP  Ham doi ngau Lagrange va bai toan doi ngau cua quy hoach toan phuong
%
%   (P)  min  0.5*x'*Q*x + c'*x   s.t.  A*x <= b,  Aeq*x = beq     (Q xac dinh duong)
%
%   L(x, lambda, mu) = 0.5*x'Q x + c'x + lambda'(A x - b) + mu'(Aeq x - beq)
%   Cuc tieu theo x:  x(lambda, mu) = -Q^{-1} r,   r = c + A'lambda + Aeq'mu
%   Ham doi ngau:     q(lambda, mu) = -0.5*r'Q^{-1} r - b'lambda - beq'mu
%   (D)  max q(lambda, mu)  s.t.  lambda >= 0
%
%   res = dual_qp(Q, c, A, b)
%   res = dual_qp(Q, c, A, b, Aeq, beq)
%
% Dau ra (struct):
%   .q       - handle q(lambda, mu)
%   .xOf     - handle x(lambda, mu), diem cuc tieu cua L theo x
%   .xP, .p  - nghiem va gia tri toi uu cua (P) (giai truc tiep)
%   .lambda, .mu, .d - nghiem va gia tri toi uu cua (D)
%   .xD      - x(lambda*, mu*): nghiem cua (P) khoi phuc tu nghiem doi ngau
%   .gap     - khe ho doi ngau p* - d*
  if nargin < 5, Aeq = []; beq = []; end
  c = c(:);
  n = numel(c);
  A = reshape(A, [], n);  b = b(:);
  Aeq = reshape(Aeq, [], n);  beq = beq(:);
  m = size(A, 1);
  p = size(Aeq, 1);
  M = [A; Aeq];
  Qi = inv(Q);

  r      = @(l, u) c + A' * l(:) + Aeq' * u(:);
  res.q  = @(l, u) -0.5 * r(l, u)' * Qi * r(l, u) - b' * l(:) - beq' * u(:);
  res.xOf = @(l, u) -Qi * r(l, u);

  % (P): giai truc tiep
  [res.xP, res.p] = solve_qp(Q, c, A, b, Aeq, beq);

  % (D): max q  <=>  min 0.5*y'G y + h'y  voi y = [lambda; mu], lambda >= 0
  G = M * Qi * M';
  h = M * Qi * c + [b; beq];
  G = (G + G') / 2;
  [y, val] = solve_qp(G, h, [], [], [], [], [zeros(m, 1); -inf(p, 1)]);
  res.lambda = y(1:m);
  res.mu     = y(m+1:end);
  res.d      = -val - 0.5 * c' * Qi * c;
  res.xD     = res.xOf(res.lambda, res.mu);
  res.gap    = res.p - res.d;
end
