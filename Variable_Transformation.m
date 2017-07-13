% Script to perform Variable Transformation.

% Written by Mayur Venkatram Lakshmi (May 2017)
% Imperial College London - Department of Aeronautics

%% Determing the transformation matrix T.
L_ij_sym = 0.5*(L_ij + L_ij');  % L_ij_sym = 0.5*(L_ij + L_ji).
A = 4*L_ij_sym;                 % A = 2*(L_ij + L_ji). 
[V,LAMBDA] = eig(A);            % Solving eigenvalue problem: A*V = LAMBDA*V.

% Eigenvalues to be sorted in descending order. lambda_1 is the largest
% eigenvalue. Find corresponing eigenvectors in correct order. The first 
% column of eigenvectors consists of the eigenvector corresponding 
% to the largest eigenvalue.
eigenvectors = fliplr(V);
disp('EIGENVECTORS of 2(L_ij + L_ji):')
disp(eigenvectors)

% Diagonalisation of A matrix: to get matrix of eigenvalues of A.
eigenvalues = eigenvectors'*A*eigenvectors;
disp('EIGENVALUES of 2(L_ij + L_ji):')
disp(eigenvalues)

% Matrices used for variable transformation.
T = eigenvectors;
T_inv = inv(eigenvectors);

%% Determining transformed matrices N_hat_xyz, N_hat_xzy, L_hat_xy and B_hat_x.
% The _hat in variable name indicates matrices in new (transformed)
% coordinates.
N_hat_xyz = zeros(N,N,N);
L_hat_xy = zeros(N,N);

% N_hat_xyz is the new N_ijk matrix after variable transformation. The nested for 
% loops implement Einstein summation over indices i,j and k.
for x = 1:N
    for y = 1:N
        for z = 1:N
            for i = 1:N
                for j = 1:N
                    for k = 1:N
                        N_hat_xyz(x,y,z) = N_hat_xyz(x,y,z) + ...
                            T_inv(x,i)*N_ijk(i,j,k)*T(j,y)*T(k,z);
                    end
                end
            end
        end
    end
end

% Similar procedure to find L_hat_xy, the L_ij matrix after coordinate 
% transformation. The nested for loops implement Einstein summation over
% indices i and j.
for x = 1:N
    for y = 1:N
        for i = 1:N
            for j = 1:N
                L_hat_xy(x,y) = L_hat_xy(x,y) + T_inv(x,i)*L_ij(i,j)*T(j,y);
            end
        end
    end
end

% B_hat_x is B_i after variable transformation. B_i is pre-multiplied by
% the inverse of the transformation matrix.
B_hat_x = T\B_i;

% Preallocation of matrix of zeros.
N_hat_xzy = zeros(N,N,N);

% Determining N_hat_xzy from N_hat_xyz.
for x = 1:N
    N_hat_xzy(x,:,:) = permute(N_hat_xyz(x,:,:),[1 3 2]);
end

%% Storing the above 3D and 2D arrays as row cell arrays.
% In the following row cell arrays, each cell corresponds to an x index value.
N_hat_xyz_array = cell(1,N);
for x = 1:N
    N_hat_xyz_array{1,x} = zeros(N,N);
    for y = 1:N
        N_hat_xyz_array{1,x}(y,:) = N_hat_xyz(x,y,:);
    end
end

N_hat_xzy_array = cell(1,N);
for x = 1:N
    N_hat_xzy_array{1,x} = zeros(N,N);
    for y = 1:N
        N_hat_xzy_array{1,x}(y,:) = N_hat_xzy(x,y,:);
    end
end

L_hat_xy_array = cell(1,N);
for x = 1:N
    L_hat_xy_array{1,x} = L_hat_xy(x,:);
end

%% Save the transformed arrays to disk.
save('transformed_arrays.mat','A','B_hat_x','B_i','L_hat_xy',...
    'L_hat_xy_array','L_ij','L_ij_sym','N','N_hat_xyz',...
    'N_hat_xyz_array','N_hat_xzy','N_hat_xzy_array','N_ijk','eigenvalues',...
    'eigenvectors','T')