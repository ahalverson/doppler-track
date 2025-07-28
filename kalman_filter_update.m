function [x_hat, P] = kalman_filter_update(x_hat_prev, P_prev, z, dt)
    % Implements a Kalman Filter update for a 2-state system.
    % State vector x = [frequency; frequency_rate]
    % Measurement z is the frequency error from the discriminator.

    % 1. State Transition Model
    % f_k = f_{k-1} + f_dot_{k-1} * dt
    % f_dot_k = f_dot_{k-1}
    A = [1 dt; 
         0 1];

    % 2. Observation Model
    % We measure the frequency directly.
    % z_k = H * x_k + v_k
    H = [1 0];

    % 3. Noise Covariances
    % Process noise: uncertainty in the motion model (e.g., satellite maneuvering)
    q_val = 0.1; % Tune this!
    Q = q_val * [dt^3/3 dt^2/2; dt^2/2 dt];
    
    % Measurement noise: uncertainty in the discriminator output
    R = 500^2; % Tune this! Depends on SNR and discriminator quality.

    % --- Prediction Step ---
    x_hat_predict = A * x_hat_prev;
    P_predict = A * P_prev * A' + Q;

    % --- Update Step ---
    % Kalman Gain
    K = P_predict * H' / (H * P_predict * H' + R);

    % Update state estimate with measurement z
    x_hat = x_hat_predict + K * (z - H * x_hat_predict);
    
    % Update the estimate covariance
    P = (eye(size(K,1)) - K * H) * P_predict;
end
