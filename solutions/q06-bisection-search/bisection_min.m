function [xmin, fmin, out] = bisection_min(f, a, b, opts)
% BISECTION_MIN  Tim cuc tieu cua ham mot bien f tren doan [a, b] bang
%                giai thuat chia doi.
%
%   [xmin, fmin, out] = bisection_min(f, a, b)
%   [xmin, fmin, out] = bisection_min(f, a, b, opts)
%
%   Hai bien the (opts.method):
%   'dichotomous' (mac dinh) - chi dung gia tri ham, f unimodal tren [a, b]:
%        x1 = c - delta, x2 = c + delta voi c = (a + b)/2;
%        f(x1) < f(x2) => cuc tieu thuoc [a, x2], nguoc lai thuoc [x1, b].
%        Moi buoc 2 lan tinh f, do dai khoang: L_{k+1} = L_k/2 + delta.
%   'derivative' - dung dau cua dao ham f' (opts.df), f kha vi, lo'i
%        (hoac f' doi dau mot lan tren [a, b]):
%        f'(c) > 0 => cuc tieu thuoc [a, c]; f'(c) < 0 => thuoc [c, b];
%        f'(c) = 0 => c la diem dung. Moi buoc 1 lan tinh f', L_{k+1} = L_k/2.
%
% Dau vao:
%   f, a, b - ham (handle) va doan tim kiem, a < b
%   opts    - (tuy chon) struct:
%             .method  : 'dichotomous' (mac dinh) | 'derivative'
%             .df      : handle dao ham f' (bat buoc voi 'derivative')
%             .tol     : do dai khoang cuoi cung (mac dinh 1e-6)
%             .delta   : nua khoang cach hai diem thu (mac dinh tol/4,
%                        phai nho hon tol/2)
%             .maxIter : so buoc lap toi da (mac dinh 200)
%             .display : in tung buoc (mac dinh false)
%
% Dau ra:
%   xmin, fmin - trung diem khoang cuoi cung va gia tri f tai do
%   out        - struct: .iterations, .nfev (so lan tinh f), .ndfev (so lan
%                tinh f'), .a, .b (day cac dau mut khoang qua tung buoc),
%                .exitflag (1: dat tol; 0: het maxIter; 2: f'(c) = 0)

  if nargin < 4, opts = struct(); end
  method  = get_opt(opts, 'method', 'dichotomous');
  tol     = get_opt(opts, 'tol', 1e-6);
  delta   = get_opt(opts, 'delta', tol / 4);
  maxIter = get_opt(opts, 'maxIter', 200);
  display = get_opt(opts, 'display', false);
  if ~(a < b)
    error('bisection_min:interval', 'Can a < b.');
  end

  useDf = strcmpi(method, 'derivative');
  if useDf
    if ~isfield(opts, 'df') || isempty(opts.df)
      error('bisection_min:df', 'Bien the ''derivative'' can opts.df = f''.');
    end
    df = opts.df;
  elseif delta >= tol / 2
    error('bisection_min:delta', ...
      'Can delta < tol/2, vi do dai khoang chi giam ve 2*delta.');
  end

  out = struct('iterations', 0, 'nfev', 0, 'ndfev', 0, 'a', a, 'b', b, 'exitflag', 0);
  if display
    fprintf('%5s %14s %14s %14s %14s\n', 'k', 'a_k', 'b_k', 'L_k', 'quyet dinh');
  end

  k = 0;
  while (b - a) > tol && k < maxIter
    c = (a + b) / 2;
    if useDf
      dc = df(c);
      out.ndfev = out.ndfev + 1;
      if dc > 0
        b = c;  how = 'f''(c) > 0';
      elseif dc < 0
        a = c;  how = 'f''(c) < 0';
      else
        a = c;  b = c;  how = 'f''(c) = 0';
        out.exitflag = 2;
      end
    else
      x1 = c - delta;
      x2 = c + delta;
      f1 = f(x1);
      f2 = f(x2);
      out.nfev = out.nfev + 2;
      if f1 < f2
        b = x2;  how = 'f(x1) < f(x2)';
      else
        a = x1;  how = 'f(x1) >= f(x2)';
      end
    end
    k = k + 1;
    out.a(end+1) = a;
    out.b(end+1) = b;
    if display
      fprintf('%5d %14.8f %14.8f %14.3e %14s\n', k, a, b, b - a, how);
    end
    if out.exitflag == 2, break; end
  end

  if out.exitflag ~= 2
    out.exitflag = double((b - a) <= tol);
  end
  xmin = (a + b) / 2;
  fmin = f(xmin);
  out.nfev = out.nfev + 1;
  out.iterations = k;
end

% =====================================================================
function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
