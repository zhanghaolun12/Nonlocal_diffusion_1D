function result=FE_solution_error_infinity_norm_time_1D(uh,accurate_function,checked_time,left,right,h_partition,basis_type,derivative_degree,Gauss_point_number)

global delta;

N_partition=(right-left)/h_partition+2*ceil(delta/h_partition);
number_of_elements=N_partition;

[M_partition,T_partition]=generate_M_T_1D(left,right,h_partition,101);

if basis_type==102
    [M_basis,T_basis]=generate_M_T_1D(left,right,h_partition,102);
elseif basis_type==101
    T_basis=T_partition;
end

%Guass quadrature's points and weights on the refenrece interval [-1,1].
[Gauss_coefficient_reference_1D,Gauss_point_reference_1D]=generate_Gauss_reference_1D(Gauss_point_number);

result=0;
%Go through all elements and accumulate the error on them.
for n=1:number_of_elements
    
    vertices=M_partition(:,T_partition(:,n));
    lower_bound=min(vertices(1),vertices(2));
    upper_bound=max(vertices(1),vertices(2));
    [Gauss_coefficient_local_1D,Gauss_point_local_1D]=generate_Gauss_local_1D(Gauss_coefficient_reference_1D,Gauss_point_reference_1D,lower_bound,upper_bound);
    uh_local=uh(T_basis(:,n));
    
    temp=max(abs(feval(accurate_function,checked_time,Gauss_point_local_1D(:))-FE_solution_1D(Gauss_point_local_1D(:),uh_local,vertices,basis_type,derivative_degree)));    
    if temp>result
        result=temp;
    end
end