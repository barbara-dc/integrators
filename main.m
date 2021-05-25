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

F_sym = [-16*s(1)+12*s(2)+16*cos(t)-13*sin(t)+u; 16*s(1)-9*s(2)-11*cos(t)+9*sin(t)+u];
dFds_sym = jacobian(F_sym, s);
dFdu_sym = jacobian(F_sym, u);

matlabFunction(F_sym, 'vars', {t,s,u}, 'file', 'dynamics');
matlabFunction(dFds_sym, 'vars', {t,s,u}, 'file', 'dFds');
matlabFunction(dFdu_sym, 'vars', {t,s,u}, 'file', 'dFdu');

%:)
[s,A,B] = expl_euler(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
[x,A,B] = expl_euler_numerical_end(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)