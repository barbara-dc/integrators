function [s,A,B] = impl_euler(t0,s,u,A,B,h,n_int)
    
    tol = 1e-5;
    for i=0:(n_int-1)
        
        k_ = [1;1];
        ti = t0+h*i;
        while norm(r(ti,s,u,k_))>tol
            k_ = k_ - nablar(ti,s,u,k_).'\r(ti,s,u,k_);
        end
        
       % [~,dFds_,dFdu_] = F(ti,s,u);
        
        dkds = -nablar(ti,s,u,k_).'\drds(ti,s,u,k_);
        dkdu = -nablar(ti,s,u,k_).'\drdu(ti,s,u,k_);
        
        dphids = eye(length(s)) + h*dkds;
        dphidu = h*dkdu;
        
        
        A=dphids*A;
        B=dphids*B+dphidu;
        
        s=s+h*k_;
    end
end