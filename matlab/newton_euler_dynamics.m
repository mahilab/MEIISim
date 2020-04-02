function [tau] = newton_euler_dynamics(a, sigma, T, mass, ms, I, f_env, m_env, q_d, q_dd, g0)

    num_links = length(mass);

    % Pad vectors for notation consistency
    a     = [0, a] + 1;
    sigma = [0, sigma];
    mass  = [0; mass];
    q_d   = [0; q_d];
    q_dd  = [0; q_dd];

    f_env = [0, f_env];
    m_env = [0, m_env];
    I     = [0, I];
    ms    = [0, ms];
    T     = [0; T];

    z = [0 0 1].';
    R         = cell(num_links+1,1);
    r         = cell(num_links+1,1);
    sigma_bar = cell(num_links+1,1);
    w_last    = cell(num_links+1,1);
    w_curr    = cell(num_links+1,1);
    w_dot     = cell(num_links+1,1);
    v_dot     = cell(num_links+1,1);
    U         = cell(num_links+1,1);
    sum_f     = cell(num_links+1,1);
    sum_m     = cell(num_links+1,1);

    f       = cell(num_links+1,1);
    m       = cell(num_links+1,1);
    tau     = sym(zeros(num_links+1,1));
    sum_k_f = cell(num_links+1,1); sum_k_f(:,1) = {[0,0,0].'};
    sum_k_m = cell(num_links+1,1); sum_k_m(:,1) = {[0,0,0].'};

    w_curr{1} = [0,0,0].';
    w_dot{1}  = [0,0,0].';
    v_dot{1}  = g0;
    U{1}      = [0,0,0;0,0,0;0,0,0];

    % Forward iterations
    for i = 2:num_links+1
        R{i}         = T{i}(1:3,1:3).'; % rotation matrix
        r{i}         = T{i}(1:3,4);     % position vector
        sigma_bar{i} = (1-sigma(i));    % not sigma

        w_last{i} = R{i}*w_curr{a(i)};                                                                  % 8.17    
        w_curr{i} = w_last{i} + sigma_bar{i}*z*q_d(i);                                                  % 8.18
        w_dot{i}  = R{i}*w_dot{a(i)} + sigma_bar{i}*(z*q_dd(i) + cross(w_last{i},z)*q_d(i));            % 8.19
        v_dot{i}  = R{i}*(v_dot{a(i)}+U{a(i)}*r{i}) + sigma(i)*(z*q_dd(i)+2*cross(w_last{i},z)*q_d(i)); % 8.20
        U{i}      = crossprodmat(w_dot{i}) + crossprodmat(w_curr{i})*crossprodmat(w_curr{i});           % 8.21
        sum_f{i}  = mass(i)*v_dot{i} + U{i}*ms{i};                                                      % 8.22
        sum_m{i}  = I{i}*w_dot{i} + cross(ms{i},v_dot{i}) + cross(w_curr{i},I{i}*w_curr{i});            % 8.23
    end

    % Backwards iterations
    for i = num_links+1:-1:2
        f{i}           = sum_f{i} + sum_k_f{i} + f_env{i};        % 8.24
        f_prev{i}      = R{i}.'*f{i};                             % 8.25
        m{i}           = sum_m{i} + sum_k_m{i} + m_env{i};        % 8.26
        tau(i)         = (sigma(i)*f{i} + sigma_bar{i}*m{i}).'*z; % 8.27

        sum_k_f{a(i)} = sum_k_f{a(i)} + f_prev{i};                          % handling summation from 8.24
        sum_k_m{a(i)} = sum_k_m{a(i)} + R{i}.'*m{i} + cross(r{i},f_prev{i}); % handling summation from 8.26
    end
    
    % Resize to not include frame 0
    tau = tau(2:(num_links+1));
end
