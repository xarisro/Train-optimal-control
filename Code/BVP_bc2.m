% The boundary conditions:
function res = BVP_bc2(ya,yb)
    global x0;
    res = [ya(1) - x0(1); ya(2) - x0(2)];
end
