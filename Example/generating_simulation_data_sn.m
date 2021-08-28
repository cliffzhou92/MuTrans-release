function [x_out, t_out, lambda_out] = generating_simulation_data_sn (F, sigma, x_0, lambda_0, par)
T_max = par. T_max;
dt = par.dt;
t_out = 0:dt:T_max;
L = length(t_out);
x = x_0;
lambda = lambda_0;
x_out= [];
lambda_out = [];

for k = 1:L
dw = sqrt(dt)*randn;
lambda = lambda+dt;
b_drift = F(x,lambda);
x = x+ b_drift*dt+ sqrt(2*sigma)*dw;
x_out = [x_out,x];
lambda_out= [lambda_out,lambda];
end

end