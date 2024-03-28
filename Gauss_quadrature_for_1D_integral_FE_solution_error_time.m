function result=Gauss_quadrature_for_1D_integral_FE_solution_error_time(uh_local,accurate_function,checked_time,vertices,Gauss_coefficient_local_1D,Gauss_point_local_1D,basis_type,derivative_degree)


Gpn=length(Gauss_coefficient_local_1D);

result=0;
for i=1:Gpn
    result=result+Gauss_coefficient_local_1D(i)*(feval(accurate_function,checked_time,Gauss_point_local_1D(i))-FE_solution_1D(Gauss_point_local_1D(i),uh_local,vertices,basis_type,derivative_degree))^2;
end