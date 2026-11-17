function [alpha, fa, ga, nfev, ngev] = wolfe_line_search(f, grad, x, d, fx, gx, alpha0, c1, c2)
% WOLFE_LINE_SEARCH  Tim do dai buoc alpha > 0 theo huong giam d thoa dieu kien
%                    Wolfe manh:
%
%       f(x + a d) <= f(x) + c1 * a * g'd                 (giam du - Armijo)
%       |grad f(x + a d)' d| <= c2 * |g'd|                (do cong)
%
%   [alpha, fa, ga, nfev, ngev] = wolfe_line_search(f, grad, x, d, fx, gx)
%   [...] = wolfe_line_search(f, grad, x, d, fx, gx, alpha0, c1, c2)
%
%   Thuat toan 3.5-3.6 cua Nocedal & Wright (Numerical Optimization): thu
%   alpha0 (mac dinh 1), mo rong (nhan doi) den khi khoang [lo, hi] chac chan
%   chua buoc thoa Wolfe, roi thu hep ("zoom") bang noi suy bac hai co bao ve.
%   Dieu kien giam du duoc noi them 10*eps*|f(x)| (sai so lam tron): gan cuc
%   tieu, luong giam c1*a*g'd nho hon sai so cua f, neu so sanh dung thi khong
%   buoc nao duoc chap nhan (Hager & Zhang 2005, "approximate Wolfe").
%
% Dau vao:
%   f, grad - ham va gradient (handle);  x, d - diem va huong (g'd < 0)
%   fx, gx  - f(x), grad f(x) da tinh san
%   alpha0  - buoc thu dau tien (mac dinh 1);  c1, c2 - mac dinh 1e-4, 0.9
%
% Dau ra:
%   alpha, fa, ga - buoc tim duoc, f va grad f tai x + alpha d
%   nfev, ngev    - so lan tinh f va grad f

  if nargin < 7 || isempty(alpha0), alpha0 = 1; end
  if nargin < 8 || isempty(c1), c1 = 1e-4; end
  if nargin < 9 || isempty(c2), c2 = 0.9; end
  dphi0 = gx' * d;
  if dphi0 >= 0
    error('wolfe_line_search:direction', 'd khong phai huong giam (g''d >= 0).');
  end
  nfev = 0;  ngev = 0;
  fx = fx + 10 * eps * abs(fx);              % noi rong theo sai so lam tron
  aPrev = 0;  fPrev = fx;  dPrev = dphi0;
  a = alpha0;
  for i = 1:60
    fa = f(x + a * d);  nfev = nfev + 1;
    if fa > fx + c1 * a * dphi0 || (i > 1 && fa >= fPrev)
      [alpha, fa, ga, nf, ng] = zoom_step(f, grad, x, d, fx, dphi0, c1, c2, aPrev, fPrev, dPrev, a, fa);
      nfev = nfev + nf;  ngev = ngev + ng;
      return;
    end
    ga = grad(x + a * d);  ngev = ngev + 1;
    da = ga' * d;
    if abs(da) <= -c2 * dphi0
      alpha = a;
      return;
    end
    if da >= 0
      [alpha, fa, ga, nf, ng] = zoom_step(f, grad, x, d, fx, dphi0, c1, c2, a, fa, da, aPrev, fPrev);
      nfev = nfev + nf;  ngev = ngev + ng;
      return;
    end
    aPrev = a;  fPrev = fa;  dPrev = da;
    a = 2 * a;                                 % mo rong buoc
  end
  alpha = a;                                   % khong tim duoc (ham khong bi chan?)
end

% =====================================================================
% Thu hep khoang giua lo (thoa Armijo, f nho hon) va hi
function [alpha, fa, ga, nf, ng] = zoom_step(f, grad, x, d, fx, dphi0, c1, c2, lo, flo, dlo, hi, fhi)
  nf = 0;  ng = 0;
  for j = 1:60
    % cuc tieu cua da thuc bac hai qua (lo, flo, dlo) va (hi, fhi)
    w = hi - lo;
    den = 2 * (fhi - flo - dlo * w);
    if den > 0
      alpha = lo - dlo * w^2 / den;
    else
      alpha = lo + w / 2;
    end
    lo2 = min(lo, hi);  hi2 = max(lo, hi);  m = 0.1 * (hi2 - lo2);
    alpha = min(max(alpha, lo2 + m), hi2 - m);     % bao ve: khong sat dau mut
    fa = f(x + alpha * d);  nf = nf + 1;
    if fa > fx + c1 * alpha * dphi0 || fa >= flo
      hi = alpha;  fhi = fa;
    else
      ga = grad(x + alpha * d);  ng = ng + 1;
      da = ga' * d;
      if abs(da) <= -c2 * dphi0
        return;
      end
      if da * (hi - lo) >= 0
        hi = lo;  fhi = flo;
      end
      lo = alpha;  flo = fa;  dlo = da;
    end
  end
  ga = grad(x + alpha * d);  ng = ng + 1;       % het so lan thu: tra ve diem cuoi
end
