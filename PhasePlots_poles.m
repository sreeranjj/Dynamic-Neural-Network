%S1 = -0.89 + 0i; S2 = -0.11 - 0i; % overdamped
S1 = -0.50 + 0.50i; S2 = -0.50 - 0.50i; % underdamped
%S1 = 0 + 0.50i; S2 = 0 - 0.50i; % undamped
%S1 = +0.01 + 0.50i; S2 = +0.01 - 0.50i; % unstable

dt = .001;
tao = 10e-3;
kv = -(S1+S2);
kp = S1*S2;

x0 = [0.5;0];
x = zeros(2,15000);
x(:,1) = x0;
i=1;

while (abs(x(1,i))>tao || abs(x(2,i))>tao) && i<15000
    x_dot = [x(2,i);-(kp*x(1,i)+kv*x(2,i))];
    x(:,i+1) = x(:,i) + x_dot*dt;
    i = i + 1;
end

plot(x(1,:),x(2,:))
hold on