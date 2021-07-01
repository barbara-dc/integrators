function [s,timegrid] = expl_euler_variable(t0,x0,u,h,Tfin)

Tabs = 1e-3;
Trel = 1e-3;


        t = t0;
        timegrid = [t0];
        s = x0;
        [F_base,~,~] = F(t,s,u);

        while (t<Tfin)
           
            s_=s+h*F_base;
            [F_,~,~] = F(t+h,s_,u);
            e = norm(F_base - F_)/(2*(Tabs + Trel*max(norm(s),norm(s_))));
            
            if e<=1
                s = s_;
                t = t+h;
                F_base = F_;
                timegrid = [timegrid, t];
            end
            
            h = h*min(2,max(0.5,0.9/sqrt(e)));
            
        end

end