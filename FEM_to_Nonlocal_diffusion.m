clc;
clear;

format short e

%近场范围
global delta;
delta=0.2;

%基函数类型
basis_type=102; 

%问题求解区域：[left,right].
left=0;
right=1;

%高斯积分点个数
Gauss_point_number=7;

%初始/结束时刻
initial_t=0;
end_t=1;
dt=0.001; %时间步长

theta=0.5;

% h_partition=1/2^5; %网格步长
% [uh,Pb,Tb]=Nonlocal_diffusion_solver_1D(left,right,h_partition,basis_type,Gauss_point_number,initial_t,end_t,dt,theta); 
% infinity_error=FE_solution_error_infinity_norm_time_1D(uh,'exact_solution',end_t,left,right,h_partition,basis_type,0,Gauss_point_number)
% L2_error=FE_solution_error_time_1D(uh,'exact_solution',end_t,left,right,h_partition,basis_type,0,Gauss_point_number)
% H1_error=FE_solution_error_time_1D(uh,'exact_solution_derivative',end_t,left,right,h_partition,basis_type,1,Gauss_point_number)
% 
% exact = exact_solution(end_t,Pb)'; %真解
% 
% plot(Pb,uh,'ro',Pb,exact,'b');
% xlabel('x'); ylabel('u'); title(strcat('the exact solution and FE solution of t =',num2str(end_t)));
% legend('有限元解','真解');
% line([left left],[min([uh;exact]) max([uh;exact])],'linestyle','--');%标记计算区域
% line([right right],[min([uh;exact]) max([uh;exact])],'linestyle','--');

j=1;
fprintf('h    infinity_error  L2_error    H1_error\n');
for i = 3 : 7
    h_partition=1/2^i;
    h(j)=h_partition;
%     h_partition = 0.2*2^(1-i);
%     h(j) = 0.2*2^(1-i);
    [uh,Pb,Tb]=Nonlocal_diffusion_solver_1D(left,right,h_partition,basis_type,Gauss_point_number,initial_t,end_t,dt,theta);
    infinity_error(j)=FE_solution_error_infinity_norm_time_1D(uh,'exact_solution',end_t,left,right,h_partition,basis_type,0,Gauss_point_number);
    L2_error(j)=FE_solution_error_time_1D(uh,'exact_solution',end_t,left,right,h_partition,basis_type,0,Gauss_point_number);
    H1_error(j)=FE_solution_error_time_1D(uh,'exact_solution_derivative',end_t,left,right,h_partition,basis_type,1,Gauss_point_number);
    fprintf('1/%d \t%7.4e \t%7.4e \t%7.4e\n',1/h_partition,infinity_error(j),L2_error(j),H1_error(j));
    j=j+1;
end
for i = 2 : length(L2_error)
    R_inf_err(i-1)=log(infinity_error(i-1)/infinity_error(i))/log(2);
    R_L2_err(i-1)=log(L2_error(i-1)/L2_error(i))/log(2);
    R_H1_err(i-1)=log(H1_error(i-1)/H1_error(i))/log(2);
end

% h=[1/2^4 1/2^5 1/2^6 1/2^7 1/2^8];
% infinity_error=[5.7707e-03 1.4874e-03 3.7749e-04 9.6513e-05 2.4398e-05];
% L2_error=[2.0053e-03 5.0168e-04 1.2688e-04 3.1965e-05 8.0556e-06];
% H1_error=[1.0369e-01 5.2073e-02 2.6420e-02 1.3256e-02 6.6638e-03];
figure
loglog(h,infinity_error,'-o',h,L2_error,'-s',h,H1_error,'-^','linewidth',2);
xlabel('log(h)'); ylabel('log(error)'); title('数值误差');
hold on
loglog([1/40 1/80],[8*10^(-5) 2*10^(-5)],'--','linewidth',2);
hold on
loglog([1/40 1/80],[2*10^(-2) 1*10^(-2)],'linewidth',2);
legend('L^{\infty}','L^2','H^1','O(h^2)','O(h)');

