function dxdt = BVP_ode8(time,x)
    global k1 k2 k3 ut T samples p12t p22t x1t_opt x2t_opt R Imin Imax
    
    index = max(ceil(time*samples/T),1);
    
    y = [x(1) - x1t_opt(index); x(2) - x2t_opt(index)];
    
    %v = -(1/R)*B'*P*y = -(1/R)*(p12*y1+p22*y2);
    v = -(1/R)*(p12t(index)*y(1) + p22t(index)*y(2));
    u = ut(index) + v;
    u = min(max(u,Imin),Imax);
    
    dxdt = [x(2); -k2*x(2)^2-k1*x(2)+k3*u];
end

