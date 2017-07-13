% Script to read in N_ijk, L_ij and B_i from file. 
% This script also verifies that SUM(N_ijk * a_i * a_j * a_k) = 0.

% Written by Mayur Venkatram Lakshmi (May 2017)
% Imperial College London - Department of Aeronautics

clear

%% User input - reads in N_ijk, L_ij and B_i from a file INPUT.mat.
load('INPUT.mat','N_ijk','L_ij','B_i')

% Examining N_ijk to determine the value of N - the size of the system.
Size_N_ijk = size(N_ijk);
N = Size_N_ijk(1);

%% Input Validation - SUMMATION TEST.
% Random vector of a values for summation test.
input_a_values = rand(N,1);

% SUM(N_ijk * a_i * a_j * a_k) should equal zero for all vectors a.

% Einstein summation over indices x,y and z.
input_sum_test = zeros(1,N);

for m = 1:N
    for i = 1:m
        for j = 1:m
            for k = 1:m
                input_sum_test(m) = input_sum_test(m) + N_ijk(i,j,k)*input_a_values(i)* ...
                    input_a_values(j)*input_a_values(k);
            end
        end
    end
end

if any(abs(input_sum_test) > 1e-10) == 1
    disp(['Warning: The input N_ijk may not '...
        'satisfy the condition N_ijk * a_i * a_j * a_k = 0'])
else
    disp('N_ijk * a_i * a_j * a_k = 0 (Input N_ijk is okay)')
end