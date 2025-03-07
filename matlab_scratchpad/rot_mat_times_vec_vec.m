function vec_rot = rot_mat_times_vec_vec(R, vec)
% Preallocate the output.
N = size(R, 1);
vec_rot = zeros(N, 3);

% Loop over each row (each individual rotation matrix and vector).
for i = 1:N
    % Reshape the i-th row of R into a 3x3 rotation matrix.
    M = reshape(R(i, :), 3, 3);
    % Multiply the rotation matrix with the i-th vector.
    vec_rot(i, :) = (M * vec(i, :)')';  
end