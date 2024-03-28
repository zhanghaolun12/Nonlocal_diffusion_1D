function [Pb,Tb]=generate_M_T_1D(left,right,h_partition,basis_type)

global delta;

h=h_partition;
left_boundary = left-delta;
right_boundary = right+delta;

if basis_type==101

   N=(right-left)/h;   %内部区域单元个数
   M=ceil(delta/h);        %边界区域单元个数/2
   Pb=zeros(1,N+2*M+1);
   Tb=zeros(2,N+2*M);

   for i=1:N+2*M+1
       Pb(1,i)=left-M*h+(i-1)*h;       
   end
   Pb(1)=left_boundary;
   Pb(end)=right_boundary;

   for i=1:N+2*M
       Tb(1,i)=i;    
       Tb(2,i)=i+1;
   end
   
elseif basis_type==102

   N=(right-left)/h;
   M=ceil(delta/h);
   dh=h/2;
   dN=N*2;
   dM=M*2;
   Pb=zeros(1,dN+2*dM+1);
   Tb=zeros(3,N+2*M);

   for i=1:dN+2*dM+1
       Pb(1,i)=left-M*h+(i-1)*dh;       
   end
   Pb(1)=left_boundary;
   Pb(2)=(Pb(1)+Pb(3))/2;
   Pb(end)=right_boundary;
   Pb(end-1)=(Pb(end-2)+Pb(end))/2;

   for i=1:N+2*M
       Tb(1,i)=2*i-1;    
       Tb(2,i)=2*i+1;
       Tb(3,i)=2*i;
   end
   
end