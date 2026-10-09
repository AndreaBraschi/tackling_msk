function y = linear(x, x_vel, x_acc, a, b)
   
y.pos = a * x + b;

% first derivative
y.vel = a * x_vel;

% second derivative
y.acc = a * x_acc;

end