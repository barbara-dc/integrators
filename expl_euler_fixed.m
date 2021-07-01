function s = expl_euler_fixed(t0,s,u,timegrid)
    for i=1:length(timegrid)-1
        ti = timegrid(i);

        [F_,~,~] = F(ti,s,u);
        
        h = timegrid(i+1) - ti;
        s=s+h*F_;
    end
end