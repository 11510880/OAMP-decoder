function [z] = clip(z, epsilon)

z(z < -epsilon) = -epsilon;
z(z > epsilon) = epsilon;
end
