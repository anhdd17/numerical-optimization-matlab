function [x, fval] = solve_nlp(prob, x0)
% SOLVE_NLP  Giai bai toan  min f(x)  s.t.  g(x) <= 0, h(x) = 0  bang
%            sqp (GNU Octave) hoac fmincon (MATLAB, Optimization Toolbox).
%            Dung de doi chieu ket qua giai tay trong cac vi du.
%
%   prob - struct: .f, .g, .h (handle; .g, .h co the bo trong)
  g = get_fun(prob, 'g');
  h = get_fun(prob, 'h');
  if exist('OCTAVE_VERSION', 'builtin') ~= 0
    % sqp cua Octave: rang buoc dang  h_eq(x) = 0  va  h_in(x) >= 0
    [x, fval] = sqp(x0(:), prob.f, h, @(x) -g(x));
  elseif exist('fmincon', 'file')
    nonlcon = @(x) deal(g(x), h(x));
    o = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'sqp');
    [x, fval] = fmincon(prob.f, x0(:), [], [], [], [], [], [], nonlcon, o);
  else
    error('solve_nlp:solver', 'Can sqp (Octave) hoac fmincon (MATLAB).');
  end
end

function fun = get_fun(prob, name)
  if isfield(prob, name) && ~isempty(prob.(name))
    fun = prob.(name);
  else
    fun = @(x) zeros(0, 1);
  end
end
