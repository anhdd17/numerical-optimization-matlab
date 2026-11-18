function w = lambert_w0(z)
% LAMBERT_W0  Nhanh chinh cua ham Lambert W voi z >= 0: nghiem w >= 0 cua
%             w * exp(w) = z (tinh tung phan tu bang phuong phap Newton).
%   Octave khong co san lambertw; MATLAB chi co trong Symbolic Math Toolbox.
  w = log1p(z);                       % xap xi dau: w*e^w = z  =>  w ~ log(1 + z)
  for k = 1:100
    ew = exp(w);
    dw = (w .* ew - z) ./ (ew .* (w + 1));
    w = w - dw;
    if all(abs(dw) <= 1e-15 * max(1, abs(w))), break; end
  end
end
