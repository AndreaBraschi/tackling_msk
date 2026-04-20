function [y] = eval_spline(cs, x, flag)

% This function evaluates the spline at the x points, given a set of
% coefficients. Furthermore, this function gives the user the opportunity
% to calulcate the 1st and 2nd time derivative of the evaluated cubic
% spline.

% NOTE --> this function assumes that the coefficients that the user
% provides as input are those of a Cubic spline. The subsequent time
% derivatives are calculated in closed form from that assumption.

%
% Inputs:
%   - cs: spline coefficients. Given that we're dealing with a Cubic
%   spline, 4 coeafficients are expected.
%   - x: vector of new data points in the input space where the user wants
%        to evaluate the spline.
%   - flag: 0, 1 or 2, depending on whether the user just wants the spline
%   evaluation (0), or also the 1st time derivative (1), or also the second
%   tim derivative (2).

% Analytical expressions:
%   1st derivative: y_dot = 0*x^3 + 3*a1*x^2 + 2*a2*x + a3
%   2nd derivative: y_dot_dot: 0*x^3 + 0*a1*x^2 + 6*a1*x + 2*a2


% evaluate spline at the new input points
y.pos = ppval(cs, x);

if flag == 1
    y.vel = cubicSpline_vel(cs, x);
elseif flag == 2
    y.vel = cubicSpline_vel(cs, x);
    y.acc = cubicSpline_acc(cs, x);
end

end



function vel = cubicSpline_vel(cs, x)
% Analytical expressions:
%   1st derivative: y_dot = 0*x^3 + 3*a1*x^2 + 2*a2*x + a3

a1 = cs.coefs(:, 1);
a2 = cs.coefs(:, 2);
a3 = cs.coefs(:, 3);

vel = 3 * a1 * x.^2 + 2 * a2 * x + a3;
end


function acc = cubicSpline_acc(cs, x)
% Analytical expressions:
%   2nd derivative: y_dot_dot: 0*x^3 + 0*a1*x^2 + 6*a1*x + 2*a2

a1 = cs.coefs(:, 1);
a2 = cs.coefs(:, 2);

acc = 6 * a1 * x + 2 * a2;
end
