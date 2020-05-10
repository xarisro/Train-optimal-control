% The boundary conditions:
function res = BVP_bc_ricatti(ya,yb)
    global S;
    res = [yb(1) - S(1,1); yb(2) - S(1,2); yb(3) - S(2,2)];
end

