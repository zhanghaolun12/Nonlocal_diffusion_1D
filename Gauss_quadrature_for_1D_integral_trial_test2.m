function result=Gauss_quadrature_for_1D_integral_trial_test2(coefficient_function_name,Gauss_coefficient_x,Gauss_point_x,vertices_x,trial_basis_type,trial_basis_index,trial_derivative_degree,Gauss_coefficient_y,Gauss_point_y,vertices_y,test_basis_type,test_basis_index,test_derivative_degree,k)

%Gpn: the number of the Gauss points of the Gauss formula we are using.
%Gauss_coefficient_local_1D,Gauss_point_local_1D:the Gauss coefficients and Gauss points on the local interval.

Gpn=length(Gauss_coefficient_y);

result=0;
for i = 1 : Gpn
    
    if k == 1
       result=result+Gauss_coefficient_x*Gauss_coefficient_y(i)*feval(coefficient_function_name,Gauss_point_x,Gauss_point_y(i))...
           *local_basis_1D(Gauss_point_x,vertices_x,trial_basis_type,trial_basis_index,trial_derivative_degree)...
           *local_basis_1D(Gauss_point_x,vertices_x,test_basis_type,test_basis_index,test_derivative_degree);
    elseif k == 2
       result=result+Gauss_coefficient_x*Gauss_coefficient_y(i)*feval(coefficient_function_name,Gauss_point_x,Gauss_point_y(i))...
           *local_basis_1D(Gauss_point_y(i),vertices_y,trial_basis_type,trial_basis_index,trial_derivative_degree)...
           *local_basis_1D(Gauss_point_x,vertices_x,test_basis_type,test_basis_index,test_derivative_degree);
    end
    
end    

