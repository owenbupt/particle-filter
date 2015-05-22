%%
% Define the related parameters
% Information for the whole
% sampling interval Ts
close all;
clear all;
clc;
Ts = 2;
% Sampling time
Times = 30;
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
R_min = 25000;
R_max = 100000;

% the detection probability
Pd = 0.8;
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
% plot(1:K, vx_target(:,1));
% plot(x_target(:, 1), y_target(:, 1));

%%
% Initialization
% there is some problems in this part, the particles are not genearated by
% random,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ao jiao de fen ge xian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for p = 1:N
    for i = 1:T
        S_kpi(1, p, i, :)= [x_target_hat(1, i) + 0 * randn; vx_target_hat(1, i) + 0 * randn; y_target_hat(1, i) + 0 * randn; vy_target_hat(1, i) + 0 * randn];
    end
    weight_kp(1, p) = 1 / N;
end

for k = 2:K
    disp(k);
    % Prediction
    for p = 1:N
        temp = zeros(4, T);
        for i = 1:T
            temp(:, i) = S_kpi(k - 1, p, i, :);
        end
        for i = 1:T
            S_kpi(k, p, i, :) = [1 Ts 0 0; 0 1 0 0; 0 0 1 Ts; 0 0 0 1]*temp(:, i) + [Ts^2 / 2, 0; Ts, 0; 0, Ts^2 / 2; 0 Ts] * [randn; randn];
        end
        
    end
    
    % Weight update
    % Problems: Only for the condition there are 4 sensors, try to find the matlab function to use a variables with a number
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for m = 1:M
        FalseArr = false_alarm(R_min, R_max, false_mean);
        [false_num, no_use] = size(FalseArr);

        Hy = detected(x_target, y_target, vx_target, vy_target, x_sen, y_sen, k, m, T, Pd, R_min, R_max, FalseArr);
        if m == 1
            Hy1 = Hy;
        elseif m == 2
            Hy2 = Hy;
        elseif m == 3
            Hy3 = Hy;
        else
            Hy4 = Hy;
        end
    end
    for p = 1:N
        for m = 1:M
            if m == 1
                Hy = Hy1;
            elseif m == 2
                Hy = Hy2;
            elseif m == 3
                Hy = Hy3;
            else
                Hy = Hy4;
            end
                    
            weight_pmk = weight_cal(S_kpi, x_sen, y_sen, p, m, k, Qr, Qtheta, Qrdot, T, Hy, false_rate_beta, Pd);
            weight_kp(k, p) = weight_kp(k, p) * weight_pmk;
        end
    end
    %%
    % normalization
    
    weight_kp(k, :) = weight_kp(k, :) ./ sum(weight_kp(k, :));
    
    % estimalation
    for i = 1:T
        x_target_hat(k, i) = sum(S_kpi(k, :, i, 1) .* weight_kp(k, :));
        vx_target_hat(k, i) = sum(S_kpi(k, :, i, 2) .* weight_kp(k, :));
        y_target_hat(k, i) = sum(S_kpi(k, :, i, 3) .* weight_kp(k, :));
        vy_target_hat(k, i) = sum(S_kpi(k, :, i, 4) .* weight_kp(k, :));
    end
    % resampling
    temp_wei = zeros(1, N);
    temp_S = zeros(1, N, T, 4);
    for i = 1:N
        temp_u = rand;
        for j = 1:N
            sum_wei = sum(weight_kp(k, 1:j));
            if sum_wei > temp_u
                temp_wei(1, i) = weight_kp(k, j);
                temp_S(1, i, :, :) = S_kpi(k, j, :, :);
                break;
            end
        end
    end
    
    weight_kp(k,:) = temp_wei(1, :);
    S_kpi(k, :, :, :) = temp_S(1, :, :, :);
    
end
for i = 1:T
    figure;
    plot(x_target(:, i), y_target(:, i), 'g*-', x_target_hat(:, i), y_target_hat(:, i), 'r+-');
    axis([-30000 30000 -30000 30000]);
    figure;
    plot(x_target(:, i) - x_target_hat(:, i), y_target(:, i) - y_target_hat(:, i), 'g*-');
end
