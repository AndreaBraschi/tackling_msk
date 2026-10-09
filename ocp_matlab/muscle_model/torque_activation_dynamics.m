function output = torque_activation_dynamics(e, a)

tau = 0.035;
output = (e-a)./tau;

end