function result=assemble_matrix_from_1D_integral2(coefficient_function_name,M_partition,T_partition,T_basis_trial,T_basis_test,number_of_trial_local_basis,number_of_test_local_basis,number_of_elements,matrix_size,Gauss_coefficient_reference_1D,Gauss_point_reference_1D,trial_basis_type,trial_derivative_degree,test_basis_type,test_derivative_degree)

global delta;
eps = 1e-8;

result=sparse(matrix_size(1),matrix_size(2));%初始化刚度矩阵

%内部积分类型
[Gauss_coefficient_reference_y,Gauss_point_reference_y]=generate_Gauss_reference_1D(8);

%Go through all elements.
%On each element, compute the volume integrals for all possible combinations of trial and test FE basis functions.
%Assemble the values of those volume integrals into the matrix.
for n=1:number_of_elements
   
    vertices_x=M_partition(:,T_partition(:,n)); %第n个单元的两个端点的全局坐标
    %单元的两个端点坐标
    lower_bound=min(vertices_x(1),vertices_x(2)); %单元坐标的下界
    upper_bound=max(vertices_x(1),vertices_x(2)); %单元坐标的上界
    
    [Gauss_coefficient_x,Gauss_point_x]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    number_of_Gauss_point_x=length(Gauss_coefficient_x); 
    
    %遍历所有外部积分点
    for i = 1 : number_of_Gauss_point_x
       left = Gauss_point_x(i)-delta;    %内部积分下限
       right = Gauss_point_x(i)+delta;   %内部积分上限
       left_element_index = get_location(left,M_partition,T_partition);
       right_element_index = get_location(right,M_partition,T_partition);
       if left_element_index == -1
           left_element_index = 1;
           left = M_partition(1);
       end
       if right_element_index == -1
           right_element_index = size(T_partition,2);
           right = M_partition(end);
       end
       
       for m = left_element_index : right_element_index
          vertices_y=M_partition(:, T_partition(:, m));
          if m == left_element_index
              [Gauss_coefficient_y,Gauss_point_y]=generate_Gauss_local_1D(Gauss_coefficient_reference_y,Gauss_point_reference_y,left,vertices_y(2));
          elseif m == right_element_index
              [Gauss_coefficient_y,Gauss_point_y]=generate_Gauss_local_1D(Gauss_coefficient_reference_y,Gauss_point_reference_y,vertices_y(1),right);
          elseif m == n  %具有奇异性的单元
              lower = vertices_x(1); %积分下限
              upper = vertices_x(2); %积分上限
              [Gauss_coefficient_reference_y,Gauss_point_reference_y]=generate_Gauss_reference_1D(8);
              [Gauss_coefficient_y,Gauss_point_y]=generate_Gauss_local_1D(Gauss_coefficient_reference_y,Gauss_point_reference_y,lower,upper);
              for alpha=1:number_of_trial_local_basis
                  for beta=1:number_of_test_local_basis
                      temp1=Gauss_quadrature_for_1D_integral_trial_test3(lower,upper,coefficient_function_name,Gauss_coefficient_x(i),Gauss_point_x(i),vertices_x,trial_basis_type,alpha,trial_derivative_degree,Gauss_coefficient_y,Gauss_point_y,vertices_y,test_basis_type,beta,test_derivative_degree,1);
                      result(T_basis_test(beta,n),T_basis_trial(alpha,n))=result(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp1;
                      
                      temp2=Gauss_quadrature_for_1D_integral_trial_test3(lower,upper,coefficient_function_name,Gauss_coefficient_x(i),Gauss_point_x(i),vertices_x,trial_basis_type,alpha,trial_derivative_degree,Gauss_coefficient_y,Gauss_point_y,vertices_y,test_basis_type,beta,test_derivative_degree,2);
                      result(T_basis_test(beta,n),T_basis_test(alpha,m))=result(T_basis_test(beta,n),T_basis_test(alpha,m))-temp2;
                  end
              end
              continue;
          else
              [Gauss_coefficient_y,Gauss_point_y]=generate_Gauss_local_1D(Gauss_coefficient_reference_y,Gauss_point_reference_y,vertices_y(1),vertices_y(2));
          end
          
          for alpha=1:number_of_trial_local_basis
              for beta=1:number_of_test_local_basis
                  temp1=Gauss_quadrature_for_1D_integral_trial_test2(coefficient_function_name,Gauss_coefficient_x(i),Gauss_point_x(i),vertices_x,trial_basis_type,alpha,trial_derivative_degree,Gauss_coefficient_y,Gauss_point_y,vertices_y,test_basis_type,beta,test_derivative_degree,1);
                  result(T_basis_test(beta,n),T_basis_trial(alpha,n))=result(T_basis_test(beta,n),T_basis_trial(alpha,n))+temp1;
                  
                  temp2=Gauss_quadrature_for_1D_integral_trial_test2(coefficient_function_name,Gauss_coefficient_x(i),Gauss_point_x(i),vertices_x,trial_basis_type,alpha,trial_derivative_degree,Gauss_coefficient_y,Gauss_point_y,vertices_y,test_basis_type,beta,test_derivative_degree,2);
                  result(T_basis_test(beta,n),T_basis_test(alpha,m))=result(T_basis_test(beta,n),T_basis_test(alpha,m))-temp2;
              end
          end
          
       end
    end

end