function [x, fval, out] = hooke_jeeves(f, x0, opts)
% HOOKE_JEEVES  Tim cuc tieu cua ham n bien f bang giai thuat tim kiem truc
%               tiep Hooke-Jeeves (pattern search), chi dung gia tri ham.
%
%   Moi vong lap gom:
%   1) Buoc do (exploratory move) quanh diem goc y: voi i = 1..n, thu
%      y + Delta*e_i, neu khong tot hon thi thu y - Delta*e_i; giu buoc nao
%      lam f giam.
%   2) Neu buoc do thanh cong (tim duoc x_new tot hon diem co so x_b):
%      buoc theo mau (pattern move)  y = x_new + (x_new - x_b), dat x_b = x_new,
%      roi lai do quanh y. Neu do quanh y khong tot hon x_b thi quay ve do quanh x_b.
%   3) Neu buoc do quanh x_b that bai: giam buoc Delta <- beta * Delta.
%   Dung khi Delta < tol.
%
%   [x, fval, out] = hooke_jeeves(f, x0)
%   [x, fval, out] = hooke_jeeves(f, x0, opts)
%
% Dau vao:
%   f    - ham muc tieu (handle), nhan vector cot n x 1
%   x0   - diem xuat phat (n x 1)
%   opts - (tuy chon) struct:
%          .delta   : buoc ban dau (mac dinh 0.5)
%          .beta    : he so giam buoc, 0 < beta < 1 (mac dinh 0.5)
%          .tol     : dung khi Delta < tol (mac dinh 1e-6)
%          .maxFev  : so lan tinh f toi da (mac dinh 1e5)
%          .display : in moi lan doi diem co so / giam buoc (mac dinh false)
%
% Dau ra:
%   x, fval - diem cuc tieu tim duoc va gia tri f
%   out     - struct: .nfev, .iterations (so lan doi diem co so),
%             .reductions (so lan giam buoc), .base (n x K: day diem co so),
%             .deltaHist, .fHist, .nfevHist (Delta, f va so lan tinh f tai moi
%             diem co so), .exitflag (1: Delta < tol,
%             0: het maxFev)

  if nargin < 3, opts = struct(); end
  delta   = get_opt(opts, 'delta', 0.5);
  beta    = get_opt(opts, 'beta', 0.5);
  tol     = get_opt(opts, 'tol', 1e-6);
  maxFev  = get_opt(opts, 'maxFev', 1e5);
  display = get_opt(opts, 'display', false);

  xb = x0(:);
  fb = f(xb);
  nfev = 1;
  n = numel(xb);
  out = struct('nfev', 0, 'iterations', 0, 'reductions', 0, 'base', xb, ...
               'deltaHist', delta, 'fHist', fb, 'nfevHist', 1, 'exitflag', 0);
  if display
    fprintf('%5s %-30s %14s %10s %s\n', 'k', 'diem co so x_b', 'f(x_b)', 'Delta', 'buoc');
    fprintf('%5d %-30s %14.8g %10.4g\n', 0, mat2str(xb', 6), fb, delta);
  end

  while delta >= tol && nfev < maxFev
    % Buoc do quanh diem co so
    [xn, fn, nfev] = explore(f, xb, fb, delta, n, nfev);
    if fn < fb
      % Thanh cong: lap buoc theo mau cho den khi khong con cai thien
      while true
        xp = xn + (xn - xb);                  % buoc theo mau
        xb = xn;  fb = fn;
        out.iterations = out.iterations + 1;
        out.base(:, end+1) = xb;
        out.deltaHist(end+1) = delta;
        out.fHist(end+1) = fb;
        out.nfevHist(end+1) = nfev;
        if display
          fprintf('%5d %-30s %14.8g %10.4g %s\n', out.iterations, mat2str(xb', 6), fb, delta, 'do + mau');
        end
        fp = f(xp);  nfev = nfev + 1;
        [xn, fn, nfev] = explore(f, xp, fp, delta, n, nfev);
        % Mau that bai neu khong giam f, hoac neu buoc do da dua ve lai chinh
        % x_b (chi lech do sai so lam tron - neu khong kiem tra se lap vo han)
        if fn >= fb || norm(xn - xb) <= delta / 2 || nfev >= maxFev
          break;                              % ve lai do quanh x_b
        end
      end
    else
      delta = beta * delta;                   % do that bai: giam buoc
      out.reductions = out.reductions + 1;
    end
  end

  out.exitflag = double(delta < tol);
  out.nfev = nfev;
  x = xb;
  fval = fb;
end

% =====================================================================
function [y, fy, nfev] = explore(f, y, fy, delta, n, nfev)
% Buoc do: lan luot thu +/- delta theo tung toa do, giu buoc lam f giam.
  for i = 1:n
    z = y;  z(i) = z(i) + delta;
    fz = f(z);  nfev = nfev + 1;
    if fz < fy
      y = z;  fy = fz;
    else
      z(i) = y(i) - delta;
      fz = f(z);  nfev = nfev + 1;
      if fz < fy
        y = z;  fy = fz;
      end
    end
  end
end

function v = get_opt(opts, name, default)
  if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
  else
    v = default;
  end
end
