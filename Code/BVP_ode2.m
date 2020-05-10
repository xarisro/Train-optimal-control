function dydt = BVP_ode2(time,y)
    global k1 k2 k3 ut T samples
    index = max(ceil(time*samples/T),1);
    u = ut(index);
    dydt = [y(2); -k2*y(2)^2-k1*y(2)+k3*u];
end
