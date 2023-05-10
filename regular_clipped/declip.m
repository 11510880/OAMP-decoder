function [zPost,vPost] = declip(y, zPri, vZpri, sigma, lambda)
%DECLIP declip using Bayes rule
% z ~ N(zPri, vZpri),  y = clip(z) + n, where n ~ N(0, sigma), lambda is
% the parameter for clip function

c1 = 1 / sqrt(pi) .* exp(-(lambda - y) .^2 / (2 .* sigma));
w1 = (lambda - zPri) / sqrt(2 .* vZpri);
c2 = 1 / sqrt(pi) .* exp(-(lambda + y) .^2 / (2 .* sigma));
w2 = (-lambda - zPri) / sqrt(2 * vZpri);
zStar = (zPri .* sigma + y .* vZpri)  / (sigma + vZpri);
c = (zPri.^2 .* sigma + y.^2 .* vZpri)  / (sigma + vZpri);
vZstar = vZpri .* sigma /  (sigma + vZpri);

c3 = sqrt(vZstar) / sqrt(pi .* vZpri ) .* exp(-(c - zStar.^2) / (2 .* vZstar));
w3 = (-lambda - zStar) / sqrt(2 .* vZstar);
w4 = (lambda - zStar) / sqrt(2 .* vZstar);

% integrals
d1 = sqrt(pi) / 2 .* erfc(w1);
d2 = sqrt(pi) / 2 .* (erfc(-inf) - erfc(w2));
d3 = sqrt(pi) / 2 .* (erfc(w3) - erfc(w4));

marginalLikelihood =  c1  .* d1 +  c2 .* d2 + c3 .* d3;

zPost = c1 .* (sqrt(vZpri/2) .* exp(-w1.^2) + zPri .* d1) + ...,
        c2 .* (-sqrt(vZpri/2) .* exp(-w2.^2) + zPri .* d2) + ...,
        c3 .* (sqrt(vZstar/2) .* (exp(-w3.^2) - exp(-w4.^2)) + zStar .* d3);
    
zPost = zPost ./ marginalLikelihood;
    

secondMoment =  c1 .* (vZpri .* (w1 .* exp(-w1.^2) + d1) + sqrt(2 .* vZpri) .* zPri .* exp(-w1.^2) + zPri.^2 .* d1) + ...,
                c2 .* (vZpri .* (-w2 .* exp(-w2.^2) + d2) - sqrt(2 .* vZpri) .* zPri .* exp(-w2.^2) + zPri.^2 .* d2) + ...,
                c3 .* (vZstar .* (w3 .* exp(-w1.^2) -w4 .* exp(-w4.^2) + d3) + sqrt(2 .* vZstar) * zStar .* (exp(-w3.^2)-exp(-w4.^2))+ zStar.^2 .* d3);
vPost = secondMoment ./ marginalLikelihood - zPost.^2;
vPost = sum(vPost) / length(vPost);
end

