function [x, positions] = generate_sparse_verctor(B,L)
% GENERATE_SPARSE_VERCTOR generate a sparse vector of with L sections with
% L >> B
% For each section, there is B values with only one non-zero element
x = zeros(B * L, 1);
positions = [];
for i = 0:L-1
    position = randi(B);
    %fprintf("non-zero postion %d \n", position);
    x(i*B+position) = sqrt(B);
    positions = [positions; i*B+position];
end
end

