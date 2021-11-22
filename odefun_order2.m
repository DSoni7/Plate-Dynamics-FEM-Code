function dydt=odefun_order2(t,z,m,k,F,idx,c)
dof=size(k,1);
idx=idx+1;
dydt(1:dof,1)=z(dof+1:2*dof);
dydt(dof+1:2*dof,1)=m*(F(:,idx)-c*z(dof+1:2*dof)-k*z(1:dof));