function F_sym = dynamics(t,in2,u1)
%DYNAMICS
%    F_SYM = DYNAMICS(T,IN2,U1)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    27-May-2021 17:55:41

s1 = in2(1,:);
s2 = in2(2,:);
t2 = cos(t);
t3 = sin(t);
t4 = s1.*1.6e+1;
F_sym = [s2.*1.2e+1+t2.*1.6e+1-t3.*1.3e+1-t4+u1;s2.*-9.0-t2.*1.1e+1+t3.*9.0+t4+u1];
