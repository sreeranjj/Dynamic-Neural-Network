dt = 0.25;
mu = [0.01;0.01];
b = [-1;-1];
P = 10;

alpha = [1.2,1,0.8]';
beta = [0.76,0.70]';
r = [ones(1,(500/dt)-1),sin(2*pi*(500:dt:1500)/250)];


x0 = [0;0];
x1 = zeros(2,floor(1500/dt)); x1(:,1:2) = [x0,x0];
x2 = zeros(2,floor(1500/dt)); x2(:,1:2) = [x0,x0];
xp = zeros(1,floor(1500/dt)); xp(1:2) = [0,0];
u1 = [0,0,0];

S = zeros(P,2);
e = zeros(P,1);
z_b = [0,0;0,0];

for k = 3:ceil(1500/dt)
    % DNU 1
    u1(2:end) = u1(1:end-1);
    u1(1) = tanh(x1(1,k));
    f1 = 1 - x1(1,k)^2;
    x1(:,k+1) = [x1(2,k);r(k) - tanh(b(1)*x1(1,k)) - tanh(b(2)*x1(2,k)*f1)];
    
    % Plant
    xp(k) = dt*(xp([k-1,k-2])*beta + u1*alpha);
    
    % DNU 2
    u2 = tanh(x2(1,k));
    f2 = 1 - x2(1,k)^2;
    x2(:,k+1) = [x2(2,k);xp(k) - tanh(b(1)*x2(1,k)) - tanh(b(2)*x2(2,k)*f2)];
    
    % Back propagation
    e(1:end-1) = e(2:end);
    e(end) = u1(1) - u2;
    
    f_prime = -2*x2(1,k);
    phi_x = [0,1;-(((sech(b(1)*x2(1,k)))^2)*b(1) + (sech(b(2)*x2(2,k)*f2)^2)*b(2)*x2(2,k)*f_prime),-((sech(b(2)*x2(2,k)*f2))^2)*b(2)*f2];
    phi_b = [0,0;-((sech(b(1)*x2(1,k)))^2)*x2(1,k),-((sech(b(1)*x2(2,k)*f2))^2)*x2(2,k)*f2];
    psi_x = [(sech(x2(1,k)))^2,0];
    psi_b = [0,0];
    
    S(1:end-1,:) = S(2:end,:);
    S(end,:) = psi_x*z_b +psi_b;
    
    z_b = phi_x*z_b + phi_b;
    
    % Weight updation
    grad = ((sum([e,e].*S,1))')/min(P,k);
    delta_b = mu.*grad;
    b = b +delta_b;
end

plot((dt:dt:1500),r);
hold on
plot((dt:dt:1500),xp);