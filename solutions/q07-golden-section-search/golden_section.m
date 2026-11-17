function [xmin, fmin, out] = golden_section(f, a, b, opts)
% GOLDEN_SECTION  Tim cuc tieu cua ham mot bien unimodal f tren doan [a, b]
%                 bang giai thuat lat cat vang.
%
%   [xmin, fmin, out] = golden_section(f, a, b)
%   [xmin, fmin, out] = golden_section(f, a, b, opts)
%
%   Voi tau = (sqrt(5) - 1)/2 ~ 0.618, hai diem thu
%       x1 = b - tau*(b - a),   x2 = a + tau*(b - a)
%   f(x1) < f(x2) => giu [a, x2], diem x1 cu tro thanh x2 moi;
%   nguoc lai     => giu [x1, b], diem x2 cu tro thanh x1 moi.
%   Vi tau^2 = 1 - tau, moi buoc chi can tinh f tai MOT diem moi, va
%   do dai khoang giam theo he so tau: L_k = tau^k * L_0.
%
% Dau vao:
%   f, a, b - ham (handle) va doan tim kiem, a < b
%   opts    - (tuy chon) struct:
%             .tol     : do dai khoang cuoi cung (mac dinh 1e-6)
%             .maxIter : so buoc lap toi da (mac dinh 500)
%             .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   xmin, fmin - trung diem khoang cuoi cung va gia tri f tai do
%   out        - struct: .iterations, .nfev (so lan tinh f), .a, .b (day cac
%                dau mut khoang), .exitflag (1: dat tol; 0: het maxIter)

  if nargin < 4, opts = struct(); end
  tol     = get_opt(opts, 'tol', 1e-6);
  maxIter = get_opt(opts, 'maxIter', 500);
  display = get_opt(opts, 'display', false);
  if ~(a < b)
    error('golden_section:interval', 'Can a < b.');
  end

  tau = (sqrt(5) - 1) / 2;
  x1 = b - tau * (b - a);
  x2 = a + tau * (b - a);
  f1 = f(x1);
  f2 = f(x2);
  out = struct('iterations', 0, 'nfev', 2, 'a', a, 'b', b, 'exitflag', 0);
  if display
    fprintf('%4s %12s %12s %12s %12s %12s %12s %11s\n', 'k', 'a_k', 'b_k', 'x1', 'x2', ...
      'f(x1)', 'f(x2)', 'L_k');
    fprintf('%4d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %11.3e\n', 0, a, b, x1, x2, f1, f2, b - a);
  end

  k = 0;
  while (b - a) > tol && k < maxIter
    if f1 < f2
      % cuc tieu thuoc [a, x2]; x1 cu -> x2 moi, tinh x1 moi
      b  = x2;
      x2 = x1;  f2 = f1;
      x1 = b - tau * (b - a);
      f1 = f(x1);
    else
      % cuc tieu thuoc [x1, b]; x2 cu -> x1 moi, tinh x2 moi
      a  = x1;
      x1 = x2;  f1 = f2;
      x2 = a + tau * (b - a);
      f2 = f(x2);
    end
    out.nfev = out.nfev + 1;
    k = k + 1;
    out.a(end+1) = a;
    out.b(end+1) = b;
    if display
      fprintf('%4d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f %11.3e\n', k, a, b, x1, x2, f1, f2, b - a);
    end
  end

  out.exitflag = double((b - a) <= tol);
  out.iterations = k;
  xmin = (a + b) / 2;
  fmin = f(xmin);
  out.nfev = out.nfev + 1;
end

% =====================================================================
function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
