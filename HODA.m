%==========================================================================
% Authors: Quoc-Viet PHAM
% Created: 11/17/2017
% Current: 04/05/2019
% E-mail:  vietpq90@gmail.com 
% Personal site: https://sites.google.com/view/vietpq90/
% 
% This function is to implement the Heuristic Offloading Decision Algorithm
% in single-server MEC systems. The original code is for multicell networks
%==========================================================================
function [ P, F, Z_l, Z_r, N_offloading, N_local ] = HODA( gArray, p_max, p_min, n0, f_l, f0, alpha, beta, q, lambda, kappa, W )
    %========================= Simulation PARAMETERS ==========================
    % hArray: channel gain matrix
    % p_max and p_min: maximum and minimum transit power of users, respectively
    % n0: background nosie
    % f_l: computational capability of mobile users
    % f0: maximum computational capability of the MEC server
    % alpha: input data size (in bits)
    % beta: number of CPU cycles required to accomplish the task
    % q: quote of MEC servers
    % lambda: weighted coefficients of latency (index 1) and energy consumption (index 2)
    % kappa: coefficient relating to the hardware architecture of the chip
    % Bs: bandwidth of a subchannel
    %==========================================================================

    % number of mobile users
    N = size(gArray, 1);    

    %===========Compute computation overhead due to local computing============
    T_l = beta./f_l;
    E_l = kappa.*beta.*(f_l).^2;
    Z_l = lambda(1)*T_l + lambda(2)*E_l;  

    % Simulation parameters
    N_HODA = q;             % four users are allowed to transmit at the same time
    rho = ones(1,N);        % preferences of different mobile users
    beta_t = lambda(1);     % weighted coefficient of latency
    beta_e = lambda(2);     % weighted coefficient of energy consumption
    zeta = 1;               % power amplifier efficiency of mobile users

    % Coefficients defined in Eq. (17)
    eta = rho.*beta_t*alpha./(W.*T_l);
    gamma = rho.*beta_e*alpha./(W.*E_l*zeta);
    tau = rho.*beta_t;

    % criterion to stop the bi-section method
    varepsilon = 1e-5;  
    phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);

    % Offloading decision and resource allocation are independently determined
    % in different SeNBs (MEC servers)
    P = zeros(1, N);        % declare the transmit power vector
    % F = zeros(N_HODA, 1);   % declare the computation resource vector
    F = zeros(N, 1);   % declare the computation resource vector
    N_local = zeros(1,N);   % the set of local computing users
    N_offloading = zeros(1,N);      % the set of offloading users
    Z_r = Inf(1, N);    % remote overhead
    
    % the set of mobile users associated with the eNB
    Nm = 1:N;
    
    % ====================== Transmit power control =======================
    for i = 1:length(Nm)
        if phi(gamma(Nm(i)), gArray(i), p_max, eta(Nm(i))) <= 0
            % line 3
            P(Nm(i)) = p_max;
        else
            % line 4
            p_s = p_min;    p_t = p_max;
            % line 13
            while abs(p_t - p_s) > varepsilon
                p_l = (p_s + p_t)/2;
                % line 8
                if phi(gamma(Nm(i)), gArray(i)/n0, p_l, eta(Nm(i))) <= 0
                    p_s = p_l;  % line 9
                else
                    p_t = p_l;  % line 11
                end
            end
            P(Nm(i)) = (p_s + p_t)/2;
        end
    end
    
    %=============== Heuristic Offloading Decision Algorithm ==============
    N_remaining_pre = Nm;
    %=========================== HODA - STAGE 1 ===========================
    for i = 1:length(Nm)
        delta_i = (eta(Nm(i)) + gamma(Nm(i))*P(Nm(i)))/log2(1 + P(Nm(i))*gArray(i)/n0) + tau(Nm(i))*f_l(Nm(i))/f0;     % Eq. (35)
        iValue = rho(Nm(i))*(beta_t + beta_e) - delta_i;   % Eq. (34)   
        if iValue <= 0  % Eq. (34) - Condition 1
            N_local(Nm(i)) = Nm(i);
            N_remaining_pre(i) = 0;
        end
    end
    
    %=========================== HODA - STAGE 2 ===========================
    % Determine S_remote and S_search - line 8 in Algorithm 2
    N_remaining = N_remaining_pre(N_remaining_pre > 0);
    N_remote = zeros(1, length(N_remaining));
    N_search_pre = N_remaining;
    for i = 1:length(N_remaining)
        gArrayIndex = (Nm == N_remaining(i));
        % Eq. (36)
        delta_i = (eta(N_remaining(i)) + gamma(N_remaining(i))*P(N_remaining(i)))/log2(1 + P(N_remaining(i))*gArray(gArrayIndex)/n0) + tau(N_remaining(i))*f_l(N_remaining(i))/f0;
        A_i_max = N_remaining;      % this set is determined after Condition 1
        A_i_max(i) = [];            % i-th should be execluded from the set
        delta_i_A = 2*sqrt(tau(N_remaining(i))*f_l(N_remaining(i)))*sum(sqrt(tau(A_i_max).*f_l(A_i_max)))/f0;
        iValue = rho(N_remaining(i))*(beta_t + beta_e) - delta_i_A - delta_i;   % Eq. (34)
        if iValue >= 0              % Eq. (34) - Condition 2
            N_remote(i) = N_remaining(i);
            N_search_pre(i) = 0;
        end
    end
    N_remote = N_remote(N_remote > 0);
    N_search = N_search_pre(N_search_pre > 0);
    % decide offloading decision for users in N_search
    N_H = N_remote; % line 9
    % if the number of users that can offload is greater than the quota
    if length(N_H) >= N_HODA    % line 10
        while length(N_H)> N_HODA   % line 11
            % line 12
            v_function = zeros(1, length(N_H));
            f_function = zeros(1, length(N_H));
            for j = 1:length(N_H)
                % allocation of the computing resources -  Eq. (30)
                f_function(j) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
                gArrayIndex = (Nm == N_H(j));
                % utility of j-th user is according to Eq. (16) and Eq. (17)
                v_function(j) = rho(N_H(j))*(beta_t + beta_e) - (eta(N_H(j)) + gamma(N_H(j))*P(N_H(j)))/log2(1 + P(N_H(j))*gArray(gArrayIndex)/n0) - tau(N_H(j))*f_l(N_H(j))/f_function(j);
            end
            % line 13
            N_H(v_function == min(v_function)) = [];
        end
         % Update the computation resources allocated by the eNB
         for j = 1:length(N_H)
             % [jk,~] = find(Nm == N_H(j));
             % F(jk) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
             F(N_H(j)) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
             gArrayIndex = (Nm == N_H(j));
             
             % update the computation overhead of offloading users
             T_t_s = alpha/(W*log2(1 + P(N_H(j))*gArray(gArrayIndex)/n0));   % Eq. (6)
             T_e_s = beta/F(N_H(j));          % Eq. (7)
             E_r = P(N_H(j))*T_t_s/zeta;    % Eq. (8)
             Z_r(N_H(j)) = lambda(1)*(T_t_s + T_e_s) + lambda(2)*E_r;
         end
    % if the number of offloading users is less than the quota - line 10
    else    % line 15
        % the while loop runs until N_H is inextensible
        iInextensible = 1;
        FLAG_exittheloop = 0;
        N_search_tmp = N_search;
        while iInextensible && ~isempty(N_search) && FLAG_exittheloop < N+1
            % line 17
            v_function = zeros(1, length(N_search));
            f_function = zeros(1, length(N_search));
            for j = 1:length(N_search)
                % allocation of the computing resources -  Eq. (30)
                f_function(j) = f0*sqrt(tau(N_search(j))*f_l(N_search(j)))/sum(sqrt(tau(N_search).*f_l(N_search)));
                gArrayIndex = (Nm == N_search(j));
                % utility of j-th user is according to Eq. (16) and Eq. (17)
                v_function(j) = rho(N_search(j))*(beta_t + beta_e) - (eta(N_search(j)) + gamma(N_search(j))*P(N_search(j)))/log2(1 + P(N_search(j))*gArray(gArrayIndex)/n0) - tau(N_search(j))*f_l(N_search(j))/f_function(j);
            end
            % find the best preferable user
            iSelected = 1;
            iCount = length(N_search);  % to void the loop
            while iSelected && iCount > 0
                iBestVal = N_search(v_function == max(v_function));
                iBestIndex = find(v_function == max(v_function));
                gArrayIndex = (Nm == N_search(iBestIndex));
                delta_iBest = (eta(N_search(iBestIndex)) + gamma(N_search(iBestIndex))*P(N_search(iBestIndex)))/log2(1 + P(N_search(iBestIndex))*gArray(gArrayIndex)/n0) - tau(N_search(iBestIndex))*f_l(N_search(iBestIndex))/f0;    % Eq. (35)
                delta_iBest_A = 2*sqrt(tau(N_search(iBestIndex))*f_l(N_search(iBestIndex)))*sum(sqrt(tau(N_H).*f_l(N_H)))/f0;     % Eq. (36)
                iBestValue = rho(N_search(iBestIndex))*(beta_t + beta_e) - delta_iBest_A - delta_iBest;   % Eq. (34)
                               
                N_Tmp = [N_H iBestVal];
                v_i_ftion = zeros(1, length(N_Tmp));
                v_ftion = zeros(1, length(N_Tmp));
                f_ftion = zeros(1, length(N_Tmp));
                for  k = 1:length(N_Tmp)
                    % utility of the new set N_Tmp
                    % allocation of the computing resources -  Eq. (30)
                    f_ftion(k) = f0*sqrt(tau(N_Tmp(k))*f_l(N_Tmp(k)))/sum(sqrt(tau(N_Tmp).*f_l(N_Tmp)));
                    gArrIndex = Nm == N_Tmp(k);
                    % utility of k-th user is according to Eq. (16) and Eq. (17)
                    v_ftion(k) = rho(N_Tmp(k))*(beta_t + beta_e) - (eta(N_Tmp(k)) + gamma(N_Tmp(k))*P(N_Tmp(k)))/log2(1 + P(N_Tmp(k))*gArray(gArrIndex)/n0) - tau(N_Tmp(k))*f_l(N_Tmp(k))/f_ftion(k);
                    % offloading decision A is indecomposable if and only
                    % if there exists no element i in A satisfying (A\{i})>A
                    % Therefore, we need to iteratively check all elements
                    % in the offloading decision set A
                    N_inDecomposable = N_Tmp;   
                    N_inDecomposable(k) = [];
                    v_ft = zeros(1, length(N_inDecomposable));
                    f_ft = zeros(1, length(N_inDecomposable));
                    for ki = 1:length(N_inDecomposable)
                        % allocation of the computing resources -  Eq. (30)
                        f_ft(ki) = f0*sqrt(tau(N_inDecomposable(ki))*f_l(N_inDecomposable(ki)))/sum(sqrt(tau(N_inDecomposable).*f_l(N_inDecomposable)));
                        gAIndex = Nm == N_inDecomposable(ki);
                        % utility of ki-th user is according to Eq. (16) and Eq. (17)
                        v_ft(ki) = rho(N_inDecomposable(ki))*(beta_t + beta_e) - (eta(N_inDecomposable(ki)) + gamma(N_inDecomposable(ki))*P(N_inDecomposable(ki)))/log2(1 + P(N_inDecomposable(ki))*gArray(gAIndex)/n0) - tau(N_inDecomposable(ki))*f_l(N_inDecomposable(ki))/f_ft(ki);
                    end
                    v_i_ftion(k) = sum(v_ft);
                end
                % check whether adding user iBest to N_H is advantagous or
                % not. In addition, check whether {N_H join iBest} is
                % decomposable or not - two conditions in line 17
                if iBestValue > 0 && sum(v_i_ftion >= sum(v_ftion)) == 0
                    iSelected = 0;  % step 17 ends
                else 
                    v_function(iBestIndex) = 0; % the current iBest is not appropriate
                    iCount = iCount - 1;        % decrease the count by 1
                    % reset the best value
                    iBestIndex = [];
                    iBestVal = [];
                end
            end
            N_search(iBestIndex) = [];      % line 18
            N_H = [N_H iBestVal];           % line 19
            
            % calculate the preference of the new set N_H (S in Algorithm 2)
            v_inextensible = zeros(1, length(N_H));
            f_inextensible = zeros(1, length(N_H));
            for j = 1:length(N_H)
                % allocation of the computing resources -  Eq. (30)
                f_inextensible(j) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
                gArrayIndex = (Nm == N_H(j));
                % utility of j-th user is according to Eq. (16) and Eq. (17)
                v_inextensible(j) = rho(N_H(j))*(beta_t + beta_e) - (eta(N_H(j)) + gamma(N_H(j))*P(N_H(j)))/log2(1 + P(N_H(j))*gArray(gArrayIndex)/n0) - tau(N_H(j))*f_l(N_H(j))/f_inextensible(j);
            end
            % to check whether S is inextensible or not
            v_inext = zeros(1, length(N_H));
            for k = 1:length(N_search)
                N_Tmp = [N_H N_search(k)];
                v_in = zeros(1, length(N_Tmp));
                f_in = zeros(1, length(N_Tmp));
                for ki = 1:length(N_Tmp)
                    % allocation of the computing resources -  Eq. (30)
                    f_in(ki) = f0*sqrt(tau(N_Tmp(ki))*f_l(N_Tmp(ki)))/sum(sqrt(tau(N_Tmp).*f_l(N_Tmp)));
                    gArrayIndex = (Nm == N_Tmp(ki));
                    % utility of ki-th user is according to Eq. (16) and Eq. (17)
                    v_in(ki) = rho(N_Tmp(ki))*(beta_t + beta_e) - (eta(N_Tmp(ki)) + gamma(N_Tmp(ki))*P(N_Tmp(ki)))/log2(1 + P(N_Tmp(ki))*gArray(gArrayIndex)/n0) - tau(N_Tmp(ki))*f_l(N_Tmp(ki))/f_in(ki);
                end
                v_inext(k) = sum(v_in);
            end
            if sum(v_inext >= sum(v_inextensible)) == 0
                iInextensible = 0;  % line 20 - S is inextensible
            end 
            
            % a trick to exit the forever loop
            if isempty(setdiff(N_search,N_search_tmp))
                FLAG_exittheloop = FLAG_exittheloop + 1;
            end
        end
        for j = 1:length(N_H)
            % [jk,~] = find(Nm == N_H(j));
            gArrayIndex = (Nm == N_H(j));
            % Update the computation resources allocated by the eNB
            % F(jk) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
            F(N_H(j)) = f0*sqrt(tau(N_H(j))*f_l(N_H(j)))/sum(sqrt(tau(N_H).*f_l(N_H)));
            % update the computation overhead of offloading users       
            T_t_s = alpha/(W*log2(1 + P(N_H(j))*gArray(gArrayIndex)/n0));   % Eq. (6)
            T_e_s = beta/F(N_H(j)); % Eq. (7)
            E_r = P(N_H(j))*T_t_s/zeta;     % Eq. (8)
            Z_r(N_H(j)) = lambda(1)*(T_t_s + T_e_s) + lambda(2)*E_r;
        end
    end
    
    % Update the set of offloading users and local users
    N_offloading(N_H) = N_H;    % offloading users
    N_local(setdiff(Nm, N_H)) = setdiff(Nm, N_H);   % local computing users

    % return the set of local of offloading users
    N_local = N_local(N_local > 0);
    N_offloading = N_offloading(N_offloading > 0);
end