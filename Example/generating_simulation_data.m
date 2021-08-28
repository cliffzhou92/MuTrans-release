function [x_out, t_out] = generating_simulation_data (F, sigma, x_0, par)

T_max = par. T_max;
dt = par.dt;
t_out = 0:dt:T_max;

L = length(t_out);
p = length(x_0);

x_out = zeros(p,L);
x_out(:,1)= x_0;

dw = sqrt(2*sigma*dt)*randn(2,L-1);

for k = 1:(L-1)
    drift = F(x_out(:,k));
    x_out(:,k+1) = x_out(:,k)+ drift*dt+ dw(:,k);
end

end