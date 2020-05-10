% The boundary conditions:
function res = BVP_bc(ya,yb)
    global x0 c1 c2 x1f;
    res = [ya(1) - x0(1); ya(2) - x0(2); yb(3) - 2*c1*(yb(1) - x1f); yb(4) - 2*c2*yb(2)];
end

