function y = cubic(x, x_vel, x_acc, a, b, c, d)
   
y.pos = a * x.^3 + b * x.^2 + c * x + d;

% first derivative
y.vel = ((3 * a * x.^2) + (2 * b * x) + c) * x_vel;

% second derivative
y.acc = ((6 * a * x) + (2 * b)) * x_vel.^2 + ...
    (((3 * a * x.^2) + (2 * b * x) + c) * x_acc);



end