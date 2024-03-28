function r=get_initial_vector(initial_function_name,M_basis)

number_of_nodes=size(M_basis,2);
r=zeros(number_of_nodes,1);
for i=1:number_of_nodes
    r(i)=feval(initial_function_name,M_basis(i));
end

end