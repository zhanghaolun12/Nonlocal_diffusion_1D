function result=assemble_vector_from_1D_integral_time(coefficient_function_name,current_time,M_partition,T_partition,T_basis_test,number_of_test_local_basis,number_of_elements,vector_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,test_basis_type,test_derivative_degree)

%vertices: the coordinates of the two vertices of a 1D element.
%Gauss_coefficient_local_1D,Gauss_point_local_1D: the Gauss coefficients and Gauss points on the local interval.


result=zeros(vector_size,1);


%Go through all elements.
%On each element, compute the 1D integrals for all test FE basis functions.
%Assemble the values of those 1D integrals into the vector.
for n=1:number_of_elements

    vertices=M_partition(:,T_partition(:,n)); %单元的两个端点坐标
    lower_bound=min(vertices(1),vertices(2)); %单元左端点坐标
    upper_bound=max(vertices(1),vertices(2)); %单元右端点坐标
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);

    for beta=1:number_of_test_local_basis     
        temp=Gauss_quadrature_for_1D_integral_test_time(coefficient_function_name,current_time,Gauss_coefficient_local_1D,Gauss_point_local_1D,vertices,test_basis_type,beta,test_derivative_degree);
        result(T_basis_test(beta,n),1)=result(T_basis_test(beta,n),1)+temp;
    end

end