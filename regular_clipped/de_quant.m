function [zPost,vPost] = de_quant(y, zPri, vZpri, sigma, quant_levels,quant_bounds)
% y = quant(z) + n
% n   N(0, sigma)
% z   N(zPri, vZpri)


L = length(quant_levels);
marginalLikelihood = 0;
zPost = 0;
secondMoment = 0;
for i=1:L
    lb = quant_bounds(i);
    rb = quant_bounds(i+1);
    quant_level = quant_levels(i);
    
    if i==1
        [res_int, coef, w] = dequant_int(quant_level, rb, y, sigma, zPri,vZpri);
%         dequant_int(quant_level, quant_bound, y, sigma, zPri,vZpri)
        % res_int = c * integ(-inf to w, exp(-t^2))
        marginalLikelihood = marginalLikelihood + res_int;
        zPost = zPost + coef .* -sqrt(vZpri/2).*exp(-w.^2) + zPri .* res_int;
        secondMoment = secondMoment + coef .* (vZpri .* -w .*exp(-w.^2) ...
            - sqrt(2.* vZpri) .* zPri .* exp(-w.^2)) ...
            + (vZpri + zPri.^2) .* res_int ;
    elseif i==L
%         [res, coef, w] = dequant_int2(quant_level, quant_bound, y, sigma, zPri,vZpri)
        [res_int, coef, w] = dequant_int2(quant_level, lb, y, sigma, zPri,vZpri);
        marginalLikelihood = marginalLikelihood + res_int;
        zPost = zPost + coef .* sqrt(vZpri/2).*exp(-w.^2) + zPri .* res_int;
        secondMoment = secondMoment + coef .* (vZpri .* w .*exp(-w.^2) ...
            + sqrt(2.* vZpri) .* zPri .* exp(-w.^2)) ...
            + (vZpri + zPri.^2) .* res_int ;
    else
        [res_int_rb, coef_rb, w_rb] = dequant_int(quant_level, rb, y, sigma, zPri,vZpri);
        [res_int_lb, coef_lb, w_lb] = dequant_int(quant_level, lb, y, sigma, zPri,vZpri);
        marginalLikelihood = marginalLikelihood + res_int_rb - res_int_lb;
        zPost = zPost + (coef_rb .* (-sqrt(vZpri/2).*exp(-w_rb.^2)) + zPri .* res_int_rb) ...
            - (coef_lb .* (-sqrt(vZpri/2).*exp(-w_lb.^2)) + zPri .* res_int_lb);
        secondMoment = secondMoment + (coef_rb .* (vZpri .* -w_rb .*exp(-w_rb.^2) ...
            - sqrt(2.* vZpri) .* zPri .* exp(-w_rb.^2)) ...
            + (vZpri + zPri.^2) .* res_int_rb) ; 
        secondMoment = secondMoment - (coef_lb .* (vZpri .* -w_lb .*exp(-w_lb.^2) ...
            - sqrt(2.* vZpri) .* zPri .* exp(-w_lb.^2)) ...
            + (vZpri + zPri.^2) .* res_int_lb) ;
    end
    
end
marginalLikelihood =   marginalLikelihood  ./ sqrt(2*pi*sigma);
zPost = zPost ./ marginalLikelihood ./ sqrt(2*pi*sigma);
secondMoment = secondMoment ./sqrt(2*pi*sigma);
vPost = secondMoment ./ marginalLikelihood - zPost.^2;
vPost = mean(vPost);

end

