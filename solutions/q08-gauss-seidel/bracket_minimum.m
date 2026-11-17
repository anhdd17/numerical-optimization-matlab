function [lo, hi, nfev] = bracket_minimum(phi, h, maxExpand)
% BRACKET_MINIMUM  Tim khoang [lo, hi] chua diem cuc tieu cua ham mot bien
%                  phi(t) quanh t = 0, bang cach mo rong buoc (nhan doi).
%
%   [lo, hi, nfev] = bracket_minimum(phi, h)
%
%   Thu t = h; neu phi tang thi doi huong (t = -h). Tiep tuc nhan doi buoc
%   cho den khi phi tang tro lai: ba diem t_{k-1} < t_k < t_{k+1} (hoac theo
%   chieu nguoc) voi phi(t_k) nho nhat cho khoang [t_{k-1}, t_{k+1}] chua
%   cuc tieu (neu phi unimodal). Neu ca hai huong deu tang, tra ve [-h, h].
  if nargin < 2 || isempty(h), h = 1; end
  if nargin < 3, maxExpand = 60; end
  f0 = phi(0);
  fp = phi(h);
  nfev = 2;
  if fp >= f0
    fm = phi(-h);
    nfev = nfev + 1;
    if fm >= f0
      lo = -h;  hi = h;                       % cuc tieu nam trong [-h, h]
      return;
    end
    h = -h;  fp = fm;                         % di theo huong am
  end
  tPrev = 0;  t = h;  ft = fp;
  for k = 1:maxExpand
    tNext = t + 2 * (t - tPrev);              % nhan doi buoc
    fNext = phi(tNext);
    nfev = nfev + 1;
    if fNext >= ft
      lo = min(tPrev, tNext);  hi = max(tPrev, tNext);
      return;
    end
    tPrev = t;  t = tNext;  ft = fNext;
  end
  error('bracket_minimum:unbounded', ...
    'Khong tim duoc khoang chua cuc tieu (ham co the khong bi chan duoi).');
end
