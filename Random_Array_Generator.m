% FORMAT: [N_ijk, L_ij, B_i, N] = Random_Array_Generator(N)

% Written by Mayur Venkatram Lakshmi (May 2017)
% Imperial College London - Department of Aeronautics

% Random_Array_Generator.m is a MATLAB function that takes as an input the
% value of N. 

% The first output is a random 3D array N_ijk (of size NxNxN) which 
% satisfies the condition N_ijk = - N_kji. The second output is a random
% 2D matrix L_ij (of size NxN). The third output is a random vector B_i of
% size Nx1, and the final output returns the value of N.

function [N_ijk, L_ij, B_i, N] = Random_Array_Generator(N)

N_ijk = zeros(N,N,N);
for j = 1:N
    A = rand(N,N);
    B = triu(A,1);
    N_ijk(:,j,:) = B  - B';
end

L_ij = rand(N,N);
B_i = rand(N,1);

end
