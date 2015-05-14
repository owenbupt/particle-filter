%%
% Define the related parameters
% Information for the whole 
% sampling interval Ts
clear all;
clc;
Ts = 150;
% Sampling time
Times = 300;
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
R_max = 10000000;

% the detection probability
Pd = 1;
%%%%%%%%%%%%%change it to 1 for debug

% false alarm probability
Pf = 0.47;
% false alarm rate
false_rate_beta = 1.0 * 10 ^ -8;
N = 2000;
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
vx = [randn, randn];
vy = [randn, randn];
Qr = 1;
Qrdot = 1;
Qtheta = 1;

for k = 1:K-1
    for i = 1:T
        x_target(k + 1, i) = x_target(k, i) + Ts * vx_target(k, i) + Ts^2 / 2 * vx(i);
        y_target(k + 1, i) = y_target(k, i) + Ts * vy_target(k, i) + Ts^2 / 2 * vy(i);
        vx_target(k + 1, i) = vx_target(k, i) + Ts * vx(i);
        vy_target(k + 1, i) = vy_target(k, i) + Ts * vy(i);
    end
end
% plot(x_target(:, 1), y_target(:, 1));

%%
% Initialization
% there is some problems in this part, the particles are not genearated by
% random,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ao jiao de fen ge xian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
weight_kp = ones(K, N);
for p = 1:N
    for i = 1:T
        S_kpi(1, p, i, :)= [x_target_hat(1, i); vx_target_hat(1, i); y_target_hat(1, i); vy_target_hat(1, i)];
    end
    weight_kp(1, p) = 1 / N;
end

for k = 2:K

% Prediction
    for p = 1:N
        temp = zeros(4, T);
        for i = 1:T
            temp(:, i) = S_kpi(k - 1, p, i, :);
        end
        for i = 1:T
            S_kpi(k, p, i, :) = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1]*temp(:, i) + [Ts^2 / 2, 0; Ts, 0; 0, Ts^2 / 2; 0 Ts] * [vx(i); vy(i)];
        end

    end

% Weight update
    for m = 1:M
        Hy = detected(x_target, y_target, vx_target, vy_target, x_sen, y_sen, k, m, T, Pd, R_min, R_max);
    end
    for p = 1:N
        for m = 1:M
            weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, M, false_rate_beta, Pd, K);
            weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
        end
    end
    %%
    % normalization
    for p = 1:N
        weight_kp(k, p) = weight_kp(k, p) / sum(weight_kp(k, :));
    end

    % resampling
    temp_wei = zeros(1, N);
    temp_S = zeros(1, N, T, 4);
    for i = 1:N
        temp_u = rand;
        sum_wei = 0;
        for j = 1:N
            sum_wei = sum_wei + weight_kp(k, j);
            if sum_wei > temp_u
                temp_wei(1, i) = weight_kp(k, j);
                temp_S(1, i, :, :) = S_kpi(k, j, :, :);
                break;
            end
        end
    end
    weight_kp(k,:) = temp_wei(1, :);
    S_kpi(k, :, :, :) = temp_S(1, :, :, :);
    for i = 1:T
        x_target_hat(k, i) = S_kpi(k, 1, i, 1);
        y_target_hat(k, i) = S_kpi(k, 1, i, 2);
        vx_target_hat(k, i) = S_kpi(k, 1, i, 3);
        vy_target_hat(k, i) = S_kpi(k, 1, i, 4);
    end
    x_target_hat(k, 1) - x_target(k, 1)
    x_target(k, 1)
end
