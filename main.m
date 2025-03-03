clear; close all;
load('DOA_Data.mat');

figure;
title('DOA with noise');
plot(theta_s_noisy, 'b-');
hold on
plot(theta_i_noisy, 'r-');
legend({'\theta_s(t)', '\theta_i(t)'}, 'Location', 'northeast');
hold off
%% 2. denoise using Kalman filter
num_steps = length(theta_s_noisy);
Q1_n = 0.01;
Q2_n = 4;


theta_s_n1 = theta_s_noisy(1);
theta_s_hat = zeros(1, 2000);
K_n1_n = 1;
F_n1_n = 1;
C_n = 4;

for n = 1 : 2000
    %G_n = Kalman_gain_computer(K_n1_n, F_n1_n, C_n, Q2_n);
    G_n = F_n1_n * K_n1_n * C_n' / (C_n * K_n1_n * C_n' + Q2_n);
    %theta_s_hat = one_step_predictor(theta_s_hat, theta_s_noisy(n), F_n1_n, C_n, G_n);
    alpha = theta_s_noisy(n) - C_n * theta_s_n1;
    theta_s_n1 = F_n1_n * theta_s_n1 + G_n * alpha;
    %K_n1_n = Riccati_equation_solver(K_n1_n, F_n1_n, C_n, G_n, Q1_n);
    K_n = K_n1_n - F_n1_n \ G_n * C_n * K_n1_n;
    K_n1_n = F_n1_n * K_n * F_n1_n' + Q1_n;
    
    theta_s_hat(n) = F_n1_n \ theta_s_n1;
end

theta_i_n1 = theta_i_noisy(1);
theta_i_hat = zeros(1, 2000);
K_n1_n = 1;
F_n1_n = 1;
C_n = 1;

for n = 1 : 2000
    %G_n = Kalman_gain_computer(K_n1_n, F_n1_n, C_n, Q2_n);
    G_n = F_n1_n * K_n1_n * C_n' / (C_n * K_n1_n * C_n' + Q2_n);
    %theta_i_hat = one_step_predictor(theta_i_hat, theta_i_noisy(n), F_n1_n, C_n, G_n);
    alpha = theta_i_noisy(n) - C_n * theta_i_n1;
    theta_i_n1 = F_n1_n * theta_i_n1 + G_n * alpha;
    %K_n1_n = Riccati_equation_solver(K_n1_n, F_n1_n, C_n, G_n, Q1_n);
    K_n = K_n1_n - F_n1_n \ G_n * C_n * K_n1_n;
    K_n1_n = F_n1_n * K_n * F_n1_n' + Q1_n;
    
    theta_i_hat(n) = F_n1_n \ theta_i_n1;
end


%% 3. plot DOAs
figure;
title('DOA denoised using Kalman filter');
plot(theta_s_hat, 'b-');
hold on
plot(theta_i_hat, 'r-');

legend({'\theta_s(t)', '\theta_i(t)'}, 'Location', 'northeast');
ylim([-20, 10]);
hold off
%% 4. beamformer using MVDR & LCMV
R = zeros(12, 12);
N = 12;
x = matX;
for n = 1 : 2000
    R = R + x(:, n) * x(:, n)';
end
R = (1 / 2000) * R + 0.01*eye(12);

theta_s_rad = deg2rad(theta_s_hat);
theta_i_rad = deg2rad(theta_i_hat);

g = [1; 1e-4];
s_hat_MVDR = zeros(1, 2000);
s_hat_LCMV = zeros(1, 2000);
w_LCMV = zeros(12, 2000);
w_MVDR = zeros(12, 2000);
a_theta_s = zeros(N, 2000);
a_theta_i = zeros(N, 2000);

for n = 1 : 2000
    a_theta_s(:, n) = exp(1j * pi * sin(theta_s_rad(n)) * (0 : N - 1)).';
    a_theta_i(:, n) = exp(1j * pi * sin(theta_i_rad(n)) * (0 : N - 1)).';
    C = [a_theta_s(:, n), a_theta_i(:, n)];
    w_MVDR(:, n) = inv(R) * a_theta_s(:, n) / (a_theta_s(:, n)' * inv(R) * a_theta_s(:, n));
    w_LCMV(:, n) = inv(R) * C * inv(C' * inv(R) * C) * g;
    s_hat_LCMV(n) = w_LCMV(:, n)' * x(:, n);
    s_hat_MVDR(n) = w_MVDR(:, n)' * x(:, n);
end


figure;
subplot(2, 1, 1);
plot(real(s_hat_MVDR));
title('Real part of $\hat{s}(t)$ for MVDR', 'Interpreter', 'latex');
subplot(2, 1, 2);
plot(imag(s_hat_MVDR));
title('Imaginary part of $\hat{s}(t)$ for MVDR', 'Interpreter', 'latex');

figure;
subplot(2, 1, 1);
plot(real(s_hat_LCMV));
title('Real part of $\hat{s}(t)$ for LCMV', 'Interpreter', 'latex');
subplot(2, 1, 2);
plot(imag(s_hat_LCMV));
title('Imaginary part of $\hat{s}(t)$ for LCMV', 'Interpreter', 'latex');

% 預先配置儲存空間
theta = linspace(-pi / 2, pi / 2, 100);
beampattern_MVDR = zeros(2000, 100);
beampattern_LCMV = zeros(2000, 100);

%% 4 - 1 MVDR beampattern

for n = 1 : 2000
    for k = 1 : 100
        %theta = theta_grid(k);
        a_theta = exp(1j * pi * sin(theta(k)) * (0 : N - 1)).';
        beampattern_MVDR(n, k) = 20 * log(abs(w_MVDR(:, n)' * a_theta));
    end
end

figure;

[Theta_display, Time] = meshgrid(theta, 1 : 2000);

surf(Theta_display, Time, beampattern_MVDR);
shading interp;
xlabel('Angle (rad)');
ylabel('Time (snapshot index)');
zlabel('Beam Pattern Magnitude (unnormalized)');
title('3D Beampattern Over Time and Angle for MVDR');
colorbar;


%% 4 - 2 LCMV beampattern

for n = 1 : 2000
    for k = 1 : 100
        a_theta = exp(1j * pi * sin(theta(k)) * (0 : N - 1)).';
        beampattern_LCMV(n, k) = 20 * log(abs(w_LCMV(:, n)' * a_theta));
    end
end


figure;
[Theta_display, Time] = meshgrid(theta, 1 : 2000);

surf(Theta_display, Time, beampattern_LCMV);
shading interp;
xlabel('Angle (rad)');
ylabel('Time (snapshot index)');
zlabel('Beam Pattern Magnitude (unnormalized)');
title('3D Beampattern Over Time and Angle for LCMV');
zlim([-100, 50]);
colorbar;


%% Robust Beamforming
% 計算資料協方差矩陣R
R = zeros(12, 12);
N = 12;
x = matX;
for n = 1 : 2000
    R = R + x(:, n) * x(:, n)';
end
R = (1 / 2000) * R;

[U, D] = eig(R);

[evl_sorted, sort_idx] = sort(diag(D), 'descend');

% 按排序重新排列特徵值和特徵向量
D_sorted = diag(evl_sorted); % 重新組成對角矩陣
U_sorted = U(:, sort_idx);          % 按列重新排列特徵向量

Q = 12;
J = 2;
e_s = U(:, 1 : J);
e_n = U(:, J + 1 : Q);

lambda = evl_sorted;

sigma_n2 = 0;
for k = J + 1 : Q
    sigma_n2 = sigma_n2 + lambda(k);
end
sigma_n2 = (1 / (Q - J)) * sigma_n2;

sigma_i2 = (lambda(2) - sigma_n2) / Q;

s_t_hat = zeros(1, 2000);
w_eig = zeros(12, 2000);

for n = 1 : 2000
    % 計算此時刻的 steering vectors
    a_theta_s = exp(1j * pi * sin(theta_s_rad(n)) * (0 : N - 1)).';
    a_theta_i = exp(1j * pi * sin(theta_i_rad(n)) * (0 : N - 1)).';
  
    
    R_tilde = sigma_i2 * a_theta_i * a_theta_i' + sigma_n2 * eye(Q);
    w_eig(:, n) = (inv(R_tilde) * a_theta_s) / (a_theta_s' * inv(R_tilde) * a_theta_s);   

    s_t_hat(n) = w_eig(:, n)' * x(:, n);
end

figure;
subplot(2, 1, 1);
plot(real(s_t_hat));
title('Real part of $\hat{s}(t)$ (Eigenspace)', 'Interpreter', 'latex');
subplot(2, 1, 2);
plot(imag(s_t_hat));
title('Imaginary part of $\hat{s}(t)$ (Eigenspace)', 'Interpreter', 'latex');

beampattern_eig = zeros(2000, k);
for n = 1 : 2000
    for k = 1 : 100
        a_theta = exp(1j * pi * sin(theta(k)) * (0 : N - 1)).';
        beampattern_eig(n, k) = 20 * log(abs(w_eig(:, n)' * a_theta));
    end
end

figure;
[Theta_display, Time] = meshgrid(theta, 1 : 2000);

surf(Theta_display, Time, beampattern_eig);
shading interp;
xlabel('Angle (rad)');
ylabel('Time (snapshot index)');
zlabel('Beam Pattern Magnitude (unnormalized)');
title('3D Beampattern Over Time and Angle for Eigenspace');
colorbar;
save('ASP_Final_results', 'theta_s_hat', 'theta_i_hat', 's_t_hat');