function dpdt = BVP_ode_ricatti(time,y)
    global x2t_opt T samples
    index = max(ceil(time*samples/T),1);
    x2 = x2t_opt(index);
    dpdt = [y(2)^2-2; y(2)*(x2/5 + 1/2) - y(1) + y(2)*y(3); 2*y(3)*(x2/5 + 1/2) - 2*y(2) + y(3)^2 - 2];
end

