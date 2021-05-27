function [s,A,B] = impl_euler(t0,s,u,A,B,h,n_int)
    
    tol = 1e-5;
    for i=0:(n_int-1)
        
        k_ = [1;1];
        ti = t0+h*i;
        while norm(r(ti,s,u,k_))>tol
            k_ = k_ - nablar(ti,s,u,k_).'\r(ti,s,u,k_);
        end
        
        %A=(eye(length(s))+h*dFds_)*A;
        %B=(eye(length(s))+h*dFds_)*B+h*dFdu_;
        
        s=s+h*k_;
    end
end