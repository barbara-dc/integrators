function [x_plus,A,B] = expl_euler_numerical_ind(t0,x,u,A,B,h,Tfin)
    [x_plus, timegrid] = expl_euler_variable(t0,x,u,h,Tfin);
    epsilon = 1e-8;
    nx = length(x);
    nu = length(u);
    A = zeros(nx);
    B = zeros(nx, nu);
    E = eye(nx);
    Eu = eye(nu);
    for i=1:nx
        x_ = x + epsilon*E(:,i);
        x_plus_ = expl_euler_fixed(t0,x_,u,timegrid);
        A(:,i) = (x_plus_ - x_plus) / epsilon;
    end
    for i=1:nu
        u_ = u + epsilon*Eu(:,i);
        x_plus_ = expl_euler_fixed(t0,x,u_,timegrid);
        B(:,i) = (x_plus_ - x_plus) / epsilon;
    end
end