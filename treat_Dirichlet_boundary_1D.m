function [A,b]=treat_Dirichlet_boundary_1D(Dirichlet_boundary_function_name,current_time,A,b,boundary_nodes,M_basis)
%Deal with Dirichlet boundary nodes.
%nbn: the total number of all the boundary nodes of FE.

nbn=size(boundary_nodes,2);

%Check all boundary nodes of FE.
for k=1:nbn

    if boundary_nodes(1,k)==-1 
        i=boundary_nodes(2,k);
        A(i,:)=0;
        A(i,i)=1;
        b(i,1)=feval(Dirichlet_boundary_function_name,current_time,M_basis(i));
    end

end