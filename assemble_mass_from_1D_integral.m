function r=assemble_mass_from_1D_integral(coefficient_function_name,M_partition,T_partition,T_basis_trial,T_basis_test,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,trial_basis_type,trial_derivative_degree,test_basis_type,test_derivative_degree)

r=sparse(matrix_size(1),matrix_size(2));

for n=1:number_of_elements
   
    vertices=M_partition(:,T_partition(:,n)); %第n个单元的两个端点的全局坐标
    %单元的两个节点坐标
    lower_bound=min(vertices(1),vertices(2)); %单元坐标的下界
    upper_bound=max(vertices(1),vertices(2)); %单元坐标的上界
    
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
   
    for alpha=1:number_of_trial_local_basis
       for beta=1:number_of_test_local_basis      
          temp=Gauss_quadrature_for_1D_integral_trial_test(coefficient_function_name,Gauss_coefficient_local_1D,Gauss_point_local_1D,vertices,trial_basis_type,alpha,trial_derivative_degree,vertices,test_basis_type,beta,test_derivative_degree); 
          r(T_basis_test(beta,n),T_basis_trial(alpha,n))=r(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp;
       end
    end

end

end