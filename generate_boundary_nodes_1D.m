function boundary_nodes=generate_boundary_nodes_1D(left,right,Pb,N_basis)


%The following boundary condition may change for different problems.
%All Dirichlet boundary nodes.
k = 1;
for i = 1 : N_basis + 1
    if Pb(i) <= left || Pb(i) >= right
        boundary_nodes(1,k)=-1;
        boundary_nodes(2,k)=i;
        k=k+1;
    end
end


