%% 4
clear all
%Constants
%syms c1 c2 x1f x0 k1 k2 k3 k4 R Imin Imax T
global c1 c2 x1f x0 k1 k2 k3 k4 R Imin Imax T
c1 = 1000;
c2 = 1000;
x1f = 10;
k1 = 0.5;
k2 = 0.1;
k3 = 1;
k4 = 10;
R = 0.3;
Imin = -2;
Imax = 2;
T = 10;
x0 = [0 0];

% State equations
syms x1 x2 p1 p2 u;
Dx1 = x2;
Dx2 = -k1*x2 -k2*x2^2 + k3*u;
% Cost function inside the integral
syms L;
L = k4 * x2 * u + R*u^2;
% Hamiltonian
syms H;
H = L + p1*Dx1 + p2*Dx2;
% Costate equations
Dp1 = -diff(H,x1);
Dp2 = -diff(H,x2);
% solve for control u
du = diff(H,u);
sol_u = solve(du, u);
% Substitute u to state equations
Dx2 = subs(Dx2, u, sol_u);
Dp2 = subs(Dp2, u, sol_u);

%% 5

% Initial guess for the solution
global samples
samples = 200;
time_lin = linspace(0,T,samples);
solinit = bvpinit(time_lin, [0 0 0.5 0.5]);

options = bvpset('Stats','on','RelTol',1e-1);

sol = bvp4c(@BVP_ode, @BVP_bc, solinit);

t = time_lin;
x = deval(sol,t);

% Calculate u(t) from x1,x2,p1,p2
global ut x1t_opt x2t_opt
ut = -(k3*x(4,:)+k4*x(2,:))/(2*R);
ut = min(max(ut,Imin),Imax);

x1t_opt = x(1,:);
x2t_opt = x(2,:); %for ricatti

%Calculate cost J
J = c1*(x(1,samples)-x1f)^2 + c2*x(2,samples)^2 + T*(k4*x(2,:)*ut' + R*(ut*ut'))/samples

%Plotting results
figure;
subplot(2,2,1)
plot(t, x(1,:));
xlabel('time')
title('Optimal x1(t)')
grid on

subplot(2,2,2)
plot(t, x(2,:), 'g');
xlabel('time')
title('Optimal x2(t)')
grid on

subplot(2,2,[3,4])
plot(t, ut, 'r');
xlabel('time')
title('Optimal saturated input u*')
grid on

%% 6
% same u input, different initial state
x0 = [0.4 0.6];

solinit2 = bvpinit(time_lin, [x0(1) x0(2)]);
sol = bvp4c(@BVP_ode2, @BVP_bc2, solinit2);
x = deval(sol,t);

J = c1*(x(1,samples)-x1f)^2 + c2*x(2,samples)^2 + T*(k4*x(2,:)*ut' + R*(ut*ut'))/samples

%Plotting results
figure;
subplot(2,2,1)
plot(t, x(1,:));
xlabel('time')
title('x1(t) with previous input')
grid on

subplot(2,2,2)
plot(t, x(2,:), 'g');
xlabel('time')
title('x2(t) with previous input')
grid on

subplot(2,2,[3,4])
plot(t, ut, 'r');
xlabel('time')
title('Previous optimal saturated input u*')
grid on

%% 7
% closed loop v(y) = -R^(-1)*B'*P(t)*y
% Need to calculate P(t) through the Ricatti equation
syms p11 p12 p22
global S B
R = 1;
Q = [2 0; 0 2];
S = [40 0; 0 40];

f = [x2; -k1*x2-k2*x2^2+k3*u];

A = [diff(f,x1) diff(f,x2)];
B = diff(f,u);

P = [p11 p12; p12 p22];

%Ricatti
DP = - P*A - A'*P + P*B*inv(R)*B'*P - Q;

solinit_ricatti = bvpinit(time_lin, [40 0 40]);
sol_ricatti = bvp4c(@BVP_ode_ricatti, @BVP_bc_ricatti, solinit_ricatti);

p = deval(sol_ricatti,t);

global p11t p12t p22t
p11t = p(1,:);
p12t = p(2,:);
p22t = p(3,:);

% v(y) = -R^(-1)*B'P(t)*y and y = x_real-x_opt
% v(y) = -(1/R)*(p12*y1+p22*y2)
%% 8
x0 = [0.4 0.6];
R = 1;

sol = bvp4c(@BVP_ode8, @BVP_bc2, solinit2);
x = deval(sol,t);

y = [x(1,:)-x1t_opt; x(2,:)-x2t_opt];
v = -(1/R)*(p12t.*y(1,:) + p22t.*y(2,:));
u_new = ut + v;
u_new = min(max(u_new,Imin),Imax);

J = c1*(x(1,samples)-x1f)^2 + c2*x(2,samples)^2 + T*(k4*x(2,:)*u_new' + R*(u_new*u_new'))/samples


%Plotting results
figure;
subplot(2,2,1)
plot(t, x(1,:));
xlabel('time')
title('x1(t) with new input')
grid on

subplot(2,2,2)
plot(t, x(2,:), 'g');
xlabel('time')
title('x2(t) with new input')
grid on

subplot(2,2,[3,4])
plot(t, u_new, 'r');
xlabel('time')
title('New optimal saturated input u*')
grid on