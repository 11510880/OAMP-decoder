function [z_vector, quant_levels,quant_bounds] = uniform_quantizer(z_vector, numLevels, quant_max)
%uniform quantizers
quant_delta = 2 * quant_max / numLevels;
quant_levels = zeros(numLevels, 1);
quant_bounds = zeros(numLevels+1, 1);
for i=1:numLevels
    quant_levels(i) = -quant_max + quant_delta / 2 + quant_delta * (i-1);
    quant_bounds(i) = -quant_max + quant_delta * (i-1);
end
quant_bounds(end) = -quant_max + quant_delta * numLevels;
for z_idx=1:length(z_vector)
    z = z_vector(z_idx);
    for i=1:numLevels
        if(z >= quant_bounds(i) && z <= quant_bounds(i+1))
            z_vector(z_idx) = quant_levels(i);
        end
    end
    if (z < quant_bounds(1))
        z_vector(z_idx) = quant_bounds(1);
    end
    if (z > quant_bounds(end))
        z_vector(z_idx) = quant_bounds(end);
    end
end
end

