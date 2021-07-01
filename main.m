commandwindow

s_init = [1;0];
u_init = 0;
t_init = 0;
t_fin = 0.1;
n_int = 100;
h = (t_fin-t_init)/n_int;
ns = length(s_init);
nu = length(u_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);
k = sym('k', [ns 1]);

F_sym = [-16*s(1)+12*s(2)+16*cos(t)-13*sin(t)+u; 16*s(1)-9*s(2)-11*cos(t)+9*sin(t)+u];
dFds_sym = jacobian(F_sym, s);
dFdu_sym = jacobian(F_sym, u);



matlabFunction(F_sym, 'vars', {t,s,u}, 'file', 'dynamics');
matlabFunction(dFds_sym, 'vars', {t,s,u}, 'file', 'dFds');
matlabFunction(dFdu_sym, 'vars', {t,s,u}, 'file', 'dFdu');

% only for implicit euler
r_sym = k - dynamics(t+h,s+h*k,u);
nablar_sym = jacobian(r_sym,k);
drds_sym = jacobian(r_sym,s);
drdu_sym = jacobian(r_sym,u);

matlabFunction(r_sym, 'vars', {t,s,u,k}, 'file', 'r');
matlabFunction(nablar_sym, 'vars', {t,s,u,k}, 'file', 'nablar');
matlabFunction(drds_sym, 'vars', {t,s,u,k}, 'file', 'drds');
matlabFunction(drdu_sym, 'vars', {t,s,u,k}, 'file', 'drdu');

% only for expl rk4
k1_sym = dynamics(t,s,u);
k2_sym = dynamics(t + 1/2 * h, s + 1/2 * h * k1_sym, u);
k3_sym = dynamics(t + 1/2 * h, s + 1/2 * h * k2_sym, u);
k4_sym = dynamics(t + h, s + h * k3_sym, u);
s_sym = s + 1/6 * h * (k1_sym +2*k2_sym + 2*k3_sym + k4_sym);
dsds_sym = jacobian(s_sym,s);
dsdu_sym = jacobian(s_sym,u);

matlabFunction(s_sym, 'vars', {t,s,u}, 'file', 's');
matlabFunction(dsds_sym, 'vars', {t,s,u}, 'file', 'dsds');
matlabFunction(dsdu_sym, 'vars', {t,s,u}, 'file', 'dsdu');

%:)
 [s_e,A_e,B_e] = expl_euler(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
 [s_e_num,A_e_num,B_e_num] = expl_euler_numerical(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
% 

[s_end,A_end,B_end] = expl_euler_numerical_end(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,t_fin)
[s_ind,A_ind,B_ind] = expl_euler_numerical_ind(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,t_fin)

% [s_im,A_im,B_im] = impl_euler(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
% [s_im_num,A_im_num,B_im_num] = impl_euler_numerical(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)

%s_im - s_im_num
 
 [s_rk,A_rk,B_rk] = expl_rk4(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
 [s_rk_num,A_rk_num,B_rk_num] = expl_rk4_numerical(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)

%s_rk - s_im
%s_rk - s_e
 
  tspan = [t_init t_fin];
[~,s_ode] = ode45(@(t,s) dynamics(t,s,0), tspan, s_init)
 


 