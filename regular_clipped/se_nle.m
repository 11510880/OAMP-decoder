function [vt] = se_nle(mmseB, tau)
%SE_LE State Evolution of NLE, generate an prediction of error for the
%LE part
vt = 1 / (1 / mmseB - 1 / tau);
end

