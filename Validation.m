% Script to validate the variable transformation and obtained bounds, and plots results.

% Written by Mayur Venkatram Lakshmi (May 2017)
% Imperial College London - Department of Aeronautics

% Important Note: 
% GAMMA and THETA in this MATLAB script represent GAMMA_hat and 
% THETA_hat, which are GAMMA and THETA after coordinate transformation.

%% Determing numerical values for ||THETA||, GAMMA, and their bounds.

% Assign values to the sdpvar variables.
assign(c_1,value(c_1))
assign(c_2,value(c_2))
assign(c_3,value(c_3))

% r is the number of experiments.
r = input('Enter value of r (number of numerical experiments): ');

% Preallocation of empty matrices and cell arrays.
Test_array = cell(r,2);
VALUE_norm_theta = zeros(r,1);
VALUE_norm_theta_bound = zeros(r,1);
Theta_results = zeros(r,1);
VALUE_Gamma = zeros(r,1);
VALUE_Gamma_bound = zeros(r,1);
VALUE_q_squared = zeros(r,1);

% c is the counter variable. r is the number of experiments, ie: the number
% of iterations.
for c = 1:r
    a_hat_values = rand(1,N);
    u_hat_values = a_hat_values(u_index);
    v_hat_values = a_hat_values(v_index);
    % Test_array holds the randomly generated a_hat (or u_hat and v_hat) values.
    % The first column holds the u_hat values and the second column holds the
    % v_hat values. The rows of Test_array correspond to iteration number.
    Test_array{c,1} = u_hat_values;
    Test_array{c,2} = v_hat_values;
    % Assigning u_hat and v_hat values to the symbolic YALMIP variables u and v. 
    assign(u,u_hat_values)
    assign(v,v_hat_values)
    
    % Caclating the numerical value of ||THETA||.
    VALUE_norm_theta(c) = sqrt(value(norm_Theta_sq));
    
    % Calculating the numerical value of P^{1/2}.
    VALUE_norm_theta_bound(c) = sqrt(value(P));
    
    % Vector Theta_results holds the value of the difference 
    % P^{1/2} - ||THETA|| for each iteration.
    Theta_results(c) = VALUE_norm_theta_bound(c) - VALUE_norm_theta(c);
    
    % Calculating numerical value of GAMMA.
    VALUE_Gamma(c) = 0;
    for x = n+1:N
        for y = n+1:N
            VALUE_Gamma(c) = VALUE_Gamma(c) + ...
                L_hat_xy_sym(x,y)*a_hat_values(x)*a_hat_values(y);
            % NOTE: L_hat_xy_sym is symmetric and diagonal from of L_hat_xy.
        end
    end
    
    VALUE_q_squared(c) = 0.5*(v_hat_values*v_hat_values');
    % VALUE_q_squared holds values of 0.5*(v_hat_i * v_hat_i).
    VALUE_Gamma_bound(c) = 0.5*eigenvalues(n+1,n+1)*VALUE_q_squared(c);
end

Gamma_over_q_sq = VALUE_Gamma ./ VALUE_q_squared;

% VALUE_Kappa is a vector. All elements are equal to the value of KAPPA.
% This vector is used to generate plots.
VALUE_Kappa = 0.5*eigenvalues(n+1,n+1)*ones(r,1);

%% Plotting the results

% Figuring out the length of the reference line (P^{1/2} = ||THETA||) 
% for ||THETA|| bound plot.
temp_y = max(VALUE_norm_theta);
temp_x = max(VALUE_norm_theta_bound);
temp = min(temp_x,temp_y);

% Plotting GAMMA/q^2 and KAPPA against iteration number.
figure
plot(1:r,Gamma_over_q_sq)
hold on
plot(1:r,VALUE_Kappa,'LineWidth',2)
xlabel('Iteration number','Interpreter','latex','FontSize',20)
title...
    ('Plot of $\hat{\Gamma}/\hat{q}^2$ and $\kappa$ against iteration number',...
    'Interpreter','latex','FontSize',20)
legend({'$\hat{\Gamma}/\hat{q}^2$','$\kappa$'},'Interpreter','latex',...
    'FontSize',20)
set(gca,'FontSize',20)

% Plotting P^{1/2} - ||THETA|| against iteration number.
figure
plot(1:r,Theta_results,'Color',[0 0.26 0.15])
title('Plot of $P^{1/2} - ||\hat{\Theta}||$ against iteration number',...
    'Interpreter','latex','FontSize',20)
xlabel('Iteration number','Interpreter','latex','FontSize',20)
ylabel('$P^{1/2} - ||\hat{\Theta}||$','Interpreter','latex','FontSize',20)
set(gca,'FontSize', 20)
 
% Plotting ||THETA|| against P^{1/2}.
figure
p = plot([0,temp*1.5],[0,temp*1.5],'r','LineWidth',2);
hold on
plot(VALUE_norm_theta_bound,VALUE_norm_theta,'*')
xlabel('$P^{1/2}$','Interpreter','latex','FontSize',20)
ylabel('$||\hat{\Theta}||$','Interpreter','latex','FontSize',20)
title('Plot of $||\hat{\Theta}||$ against $P^{1/2}$','Interpreter',...
    'latex','FontSize',20)
legend({'$P^{1/2} = ||\hat{\Theta}||$','Numerical Data'},'Interpreter',...
    'latex','FontSize',20)
set(gca,'FontSize',20)
grid on

% Determining maximum and minimum error in the ||THETA|| and GAMMA bounds
% over all iterations.
Gamma_max_error = max(VALUE_Kappa - Gamma_over_q_sq);
Gamma_min_error = min(VALUE_Kappa - Gamma_over_q_sq);
Theta_max_error = max(Theta_results);
Theta_min_error = min(Theta_results);

%% Determining da_i_dt and da_i_hat_dt.
% Generating random vector of a_i values. From this vector, the corresponding
% transformed a_hat vector is calculated by premultiplying by A^{-1}.
a_values = rand(N,1);
a_hat_values = T\a_values;

% Preallocation of matrix of zeros.
Term_A1 = zeros(N,1);
Term_A2 = zeros(N,1);
Term_B1 = zeros(N,1);
Term_B2 = zeros(N,1);

% Einstein summation over indices j and k.
for i = 1:N
    for j = 1:N
        for k = 1:N
            Term_A1(i) = Term_A1(i) + N_ijk(i,j,k)*a_values(j)*a_values(k);
        end
    end
end

% Summing over index j.
for i = 1:N
    for j = 1:N
        Term_A2(i) = Term_A2(i) + L_ij(i,j)*a_values(j);
    end
end

% Summing over indices y and z.
for x = 1:N
    for y = 1:N
        for z = 1:N
            Term_B1(x) = Term_B1(x) + N_hat_xyz(x,y,z)*a_hat_values(y)*a_hat_values(z);
        end
    end
end

% Summing over index y.
for x = 1:N
    for y = 1:N
        Term_B2(x) = Term_B2(x) + L_hat_xy(x,y)*a_hat_values(y);
    end
end

% If variable transformation is correct: da_i_dt = T * da_i_hat_dt.
da_i_dt = Term_A1 + Term_A2;
da_i_hat_dt = Term_B1 + Term_B2;

%% Checking if N_hat_xyz * a_hat_x * a_hat_y * a_hat_z = 0
sum_test = zeros(1,N);
for m = 1:N
    for x = 1:m
        for y = 1:m
            for z = 1:m
                sum_test(m) = sum_test(m) + N_hat_xyz(x,y,z)*a_hat_values(x)* ...
                    a_hat_values(y)*a_hat_values(z);
            end
        end
    end
end

if any(abs(sum_test) > 1e-10) == 1
    disp(['Warning: The input N_hat_xyz may not '...
        'satisfy the condition N_hat_xyz * a_hat_x * a_hat_y * a_hat_z = 0'])
else
    disp('N_hat_xyz * a_hat_x * a_hat_y * a_hat_z = 0 -- All okay')
end

%% Saving validation data in a MATLAB workspace file.
save('Validation_Data.mat','Test_array','VALUE_norm_theta','VALUE_norm_theta_bound',...
    'Theta_results','VALUE_Gamma','VALUE_q_squared','VALUE_Gamma_bound',...
    'Gamma_over_q_sq','VALUE_Kappa','Gamma_max_error','Gamma_min_error',...
    'Theta_max_error','Theta_min_error')