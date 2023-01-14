dt = .001;
tao = 10e-4;

kv = 0.5;
kp = 1;

x0 = [0.5;0]; % default initial condition

%x0 = [2.0001;0]; % limit cycle for squared controller & (kp,kv) = (1,0.5)
%x0 = [2.35;0]; % limit cycle for mod controller & (kp,kv) = (1,0.5)
%x0 = [2.006717;0]; % limit cycle for squared controller & (kp,kv) = (1,1)

x = zeros(2,50000);
x(:,1) = x0;
i=1;

while (abs(x(1,i))>tao || abs(x(2,i))>tao) && i<50000
    kv1 = kv; % linear controller
    %kv1 = kv*(1 - abs(x(1,i))); % mod controller
    %kv1 = kv*(1 - x(1,i)^2); % squared controller
    x_dot = [x(2,i);-(kp*x(1,i)+kv1*x(2,i))];
    x(:,i+1) = x(:,i) + x_dot*dt;
    i = i + 1;
end

plot(x(1,:),x(2,:))
hold on