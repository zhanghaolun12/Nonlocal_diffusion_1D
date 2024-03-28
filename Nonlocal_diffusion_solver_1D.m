function [r,M_basis,T_basis]=Nonlocal_diffusion_solver_1D(left,right,h_partition,basis_type,Gauss_point_number,initial_t,end_t,dt,theta)
% 一维非局部波动方程求解器
global delta;

N_t=(end_t-initial_t)/dt;

% left_boundary = left - delta;        %定义新的端点
% right_boundary = right + delta;
N_partition=(right-left)/h_partition+2*ceil(delta/h_partition);  %总的区间段数

if basis_type==102
    N_basis=N_partition*2;
elseif basis_type==101
    N_basis=N_partition;
end

%Mesh information for partition and finite element basis functions.
[M_partition,T_partition]=generate_M_T_1D(left,right,h_partition,101);

%FE information
if basis_type==102
    [M_basis,T_basis]=generate_M_T_1D(left,right,h_partition,102);
elseif basis_type==101
    M_basis=M_partition;
    T_basis=T_partition;
end 


%Guass quadrature's points and weights on the refenrece interval [-1,1].
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number);

%Assemble the stiffness matrix.
number_of_elements=N_partition;
matrix_size=[N_basis+1 N_basis+1]; %矩阵大小
if basis_type==102
    number_of_trial_local_basis=3;
    number_of_test_local_basis=3;
elseif basis_type==101
    number_of_trial_local_basis=2;
    number_of_test_local_basis=2;
end

%Assemble the stiffness matrix
A=assemble_matrix_from_1D_integral('function_a',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,basis_type,0);

%Assemble the mass matrix
M=assemble_mass_from_1D_integral('function_one',M_partition,T_partition,T_basis,T_basis,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0,basis_type,0);

A_tilde = M/dt+theta*A;

%Get the information matrix for boundary nodes.
boundary_nodes=generate_boundary_nodes_1D(left,right,M_basis,N_basis);

% Get the initial solution.
X_old = get_initial_vector('function_initial',M_basis);

vector_size=N_basis+1; %向量大小
b2=assemble_vector_from_1D_integral_time('function_f',initial_t,M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0);
for n = 1 : N_t
    current_time=initial_t+dt*n;
    
    %Assemble the load vector.
    b1=assemble_vector_from_1D_integral_time('function_f',current_time,M_partition,T_partition,T_basis,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,basis_type,0);
    b_tilde=theta*b1+(1-theta)*b2+(M/dt)*X_old-(1-theta)*A*X_old;

    %Deal with Dirichlet boundary condition.
    [A_tilde,b_tilde]=treat_Dirichlet_boundary_1D('function_g',current_time,A_tilde,b_tilde,boundary_nodes,M_basis);
    
    %Compute the solution.
    X=A_tilde\b_tilde;
    
    X_old=X;
    
    r=X_old;
    b2=b1;
    
    if n == N_t
        maxerror=get_maximum_error_time_1D(end_t,X_old,N_basis,left,h_partition);
    end
    
end
