function maxerror=get_maximum_error_time_1D(checked_time,solution,N_basis,left,h_basis)

maxerror=0;
for i=1:N_basis+1
    temp=solution(i)-exact_solution(checked_time,left+(i-1)*h_basis);
    if abs(maxerror)<abs(temp)
        maxerror=temp;
    end
end