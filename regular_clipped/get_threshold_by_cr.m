function [epsilon] = get_threshold_by_cr(CR,z)
%GET_THRESHOLD_BY_CR �˴���ʾ�йش˺�����ժҪ
epsilon = 10 ^ (CR / 10) * (mean(z.^2));
epsilon = 10 ^ (CR / 20);
end

