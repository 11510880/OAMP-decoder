function [res, coef, w] = dequant_int2(quant_level, quant_bound, y, sigma, zPri,vZpri)
%compute the integral for dequantizer
% res = int from bound to inf
%   此处显示详细说明
coef = 1 / sqrt(pi) .* exp(-(y-quant_level) .^2 / (2 .* sigma));
w = (quant_bound - zPri) / sqrt(2 * vZpri);
d = sqrt(pi) / 2 .* (erfc(w));
res = d .* coef;
end

