function [tau] = se_le(vt, mmseA)
%SE_LE State Evolution of LMMSE, generate an prediction of error for the
%NLE part
tau =  1 / (1 / mmseA - 1/vt);
end

