function dydt = BVP_ode(t,y)
    global k1 k2 k3 k4 R Imin Imax
    u = -(k3*y(4) + k4*y(2))/(2*R);
    u = min(max(u,Imin),Imax);
    dydt = [y(2); -k2*y(2)^2-k1*y(2)+k3*u; 0; y(4)*(k1+2*k2*y(2)) - k4*u - y(3)];
end

