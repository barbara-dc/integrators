commandwindow

s_init = [0;0;0;0];
u_init = 1;
t_init = 0;
t_fin = 0.1;
n_int = 10;
h = (t_fin-t_init)/n_int;
ns = length(s_init);
nu = length(u_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);
k = sym('k', [ns 1]);

F_sym = [s(3); 
    s(4); 
    (m*l*sin(s(2))*s(4)^2 + m*gravity*cos(s(2))*sin(s(2)) + u) / (M+m-m*(cos(s(2)))^2);
    -((m*l*cos(s(2))*sin(s(2))*s(4)^2 + (M+m)*gravity*sin(s(2)) + u*cos(s(2))) / (l*(M+m-m*(cos(s(2)))^2)))];
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
[s_eul,A_eul,B_eul] = expl_euler(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
[s_eul_num,A_eul_num,B_eul_num] = expl_euler_numerical_end(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
% 
 [s_,A_,B_] = expl_rk4(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
 [s,A,B] = expl_rk4_numerical(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
 tspan = [t_init t_fin];
[t,s_ode] = ode45(@(t,s) dynamics(t,s,1), tspan, s_init);