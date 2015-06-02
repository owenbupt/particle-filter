clear('all');
clc;

 %%
% Define the related parameters
% Information for the whole
% sampling interval Ts
Ts = 2;
% Sampling time
Times = 6;
% Sampling numbers
K = Times / Ts;
T = 2; % targets numbers
M = 4; % sensors numbers

% Information for sensors, in this example there are 4 sensors.
x_sen0 = 0;
y_sen0 = 0;

sigma_r_sen0 = 500;
sigma_theta_sen0 = 2;
sigma_rdot_sen0 = 6;

x_sen1 = -30000;
y_sen1 = 0;

sigma_r_sen1 = 500;
sigma_theta_sen1 = 2;
sigma_rdot_sen1 = 6;

x_sen2 = 15000;
y_sen2 = 26000;

sigma_r_sen2 = 500;
sigma_theta_sen2 = 2;
sigma_rdot_sen2 = 6;

x_sen3 = 15000;
y_sen3 = -26000;

sigma_r_sen3 = 250;
sigma_theta_sen3 = 1;
sigma_rdot_sen3 = 3;


x_sen = [x_sen0, x_sen1, x_sen2, x_sen3];
y_sen = [y_sen0, y_sen1, y_sen2, y_sen3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sensors sensible area
R_min = 0;
R_max = 100000000;

% the detection probability
Pd = 1;
%%%%%%%%%%%%%change it to 1 for debug

% false alarm probability
Pf = 0.47;
% false alarm rate
false_rate_beta = 1.0 * 10 ^ -8;
N = 200;
% numbers of particles


% targets Information
x_target0 = 2500;
y_target0 = 25000;
vx_target0 = 0;
vy_target0 = -222;
x_target1 = -5000;
y_target1 = -20000;
vx_target1 = 120;
vy_target1 = 0;
x_target = zeros(K, T);
vx_target = zeros(K, T);
y_target = zeros(K, T);
vy_target = zeros(K, T);
x_target(1, :) = [x_target0, x_target1];
y_target(1, :) = [y_target0, y_target1];
vx_target(1, :) = [vx_target0, vx_target1];
vy_target(1, :) = [vy_target0, vy_target1];
x_target_hat = x_target;
y_target_hat = y_target;
vx_target_hat = vx_target;
vy_target_hat = vy_target;

S_kpi = zeros(K, N, T, 4);
% vx = [randn, randn];
% vy = [randn, randn];
Qr = 1;
Qrdot = 1;
Qtheta = 1;
false_mean = 1;
weight_kp = ones(K, N);
% caculate the real location and v
for k = 1:K-1
    for i = 1:T
        a = randn;
        b = randn;
        x_target(k + 1, i) = x_target(k, i) + Ts * vx_target(k, i) + Ts^2 / 2 * a;
        vx_target(k + 1, i) = vx_target(k, i) + Ts * a;
        y_target(k + 1, i) = y_target(k, i) + Ts * vy_target(k, i) + Ts^2 / 2 * b;
        vy_target(k + 1, i) = vy_target(k, i) + Ts *  b;
    end
end

for p = 1:N
    for i = 1:T
        S_kpi(1, p, i, :)= [x_target_hat(1, i) + 0 * randn; vx_target_hat(1, i) + 0 * randn; y_target_hat(1, i) + 0 * randn; vy_target_hat(1, i) + 0 * randn];
    end
    weight_kp(1, p) = 1 / N;
end


% A = zeros(10, 4);
Achose = zeros(1, 4);
A_km = zeros(K, M);
chosn_m = zeros(1, T);
Cost_record = zeros(1, 16^K);
cost_index = 1;
S_kpi_m = zeros(16^K, N, T, 4);
weight_Kkp = zeros(16^K, K, N);
% for i = 1:10
%     a = dec2bin(i, 4);
%     A(i,:) = a;
% end
temp_S_kpi = S_kpi;
for k = 2:K
    for t = 1:16^(k - 2)
        for m1 = 1:M
            for m2 = 1:M
                A_km(k, :) = 0;
                A_km(k, m1) = 1;
                A_km(k, m2) = 1;
                if k ~= 2
                    S_kpi(k-1,:,:,:) = S_kpi_m((t + 16^(k - 3)),:,:,:);
                    weight_kp(k, :) = 1;
                    [x_target_hat, vx_target_hat, y_target_hat, vy_target_hat, S_kpi, weight_kp] = fpf(k, A_km, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat);
                    Cost_record(1, cost_index) = Cost_record(1, t + floor(16^(k - 3))) + one_step_cost(k, A_km, x_target_hat, y_target_hat, x_target, y_target);
                    
                else
                    S_kpi = temp_S_kpi;
                    weight_kp(k, :) = 1;
                    [x_target_hat, vx_target_hat, y_target_hat, vy_target_hat, S_kpi, weight_kp] = fpf(k, A_km, S_kpi, weight_kp, x_target, vx_target, y_target, vy_target, x_target_hat, vx_target_hat, y_target_hat, vy_target_hat);
                    Cost_record(1, cost_index) = one_step_cost(k, A_km, x_target_hat, y_target_hat, x_target, y_target);
                end
%                 weight_kp(k,:) = weight_Kkp(cost_index, k, :);
                S_kpi_m(cost_index,:,:,:) = S_kpi(k,:,:,:);
                cost_index = cost_index + 1;
                
            end
        end
    end

end
total_cost = min(Cost_record(1, cost_index - 16^(K-1):cost_index - 1));
disp(total_cost);

% for i = 1:T
%     figure;
%     plot(x_target(:, i), y_target(:, i), 'g*-', x_target_hat(:, i), y_target_hat(:, i), 'r+-');
%     axis([-30000 30000 -30000 30000]);
%     figure;
%     plot(1:K, x_target(:, i) - x_target_hat(:, i), 1:K, y_target(:, i) - y_target_hat(:, i));
% end



