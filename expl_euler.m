function [s,A,B] = expl_euler(t0,s,u,A,B,h,n_int)
    for i=0:(n_int-1)
        ti = t0+h*i;
        % get F,dFds,dFd
        [F_,dFds_,dFdu_] = F(ti,s,u);
        A=(eye(length(s))+h*dFds_)*A;
        B=(eye(length(s))+h*dFds_)*B+h*dFdu_;
        s=s+h*F_;
    end
end