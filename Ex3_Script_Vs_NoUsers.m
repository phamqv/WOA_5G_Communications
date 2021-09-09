% =========================================================================
% =============================== EXAMPLE 3 ===============================
% =========================================================================

% Load details of the selected benchmark function for the WOA algorithm
Function_name = 'Ex3';
% Simulation parameters
NoUsers = 2:1:10;       % number of mobile users
S = 4;                  % number of subchannels in each cell
cellRadiusMax = 250;    % maximum radius of a small cell
cellRadiusMin = 5;      % minimum radius of a small cell
logNormalMean = 0;      % 
logNormalDeviation = 8; % 
noRealization = 20;   	% number of realizations
lambda_t = 0.5;         % weighted parameter of computation time
lambda_e = 1-lambda_t;  % weighted parameter of energy consumption
lambda = [lambda_t lambda_e];
n0 = db2lin(-100 - 30); % power spectral density of additive Gaussian noise
W = 1e6;               % bandwidth of a subchannel
p_max = 100e-3;         % maximum transmit power of each mobile user
p_min = 1e-6;           % minimum transmit power of each mobile user
f0 = 4.0*1e9;          	% computational capability of the MEC server
alpha = 1*420e3;        % input data size of computation tasks = 420KB
beta = 1000e6;          % total required number of CPU cycles of mobile sers = 1000 Megacycles       
kappa = 5e-27;			% coefficient depdending on the chip's hardware architecture
zeta = 1;               % power amplifier efficiency of mobile users
% local computing capabilities
% f_l = 1e9*[0.5, 0.8, 0.8, 1, 1, 0.5];
f_local = 1e9*[0.5 0.8 1.0];
f_user = zeros(1,1000); 
for n = 1:1000
    f_user(n) = f_local(randi(length(f_local),1));
end
% parameters for the WOA algorithm
SearchAgents_no = 30;   % Number of search agents
Max_iter = 100;         % Maximum numbef of iterations
BWOA_num = 9;

% Define and declare performance metrics
swco_HODA = zeros(length(NoUsers), noRealization);
po_HODA = zeros(length(NoUsers), noRealization);
su_HODA = zeros(length(NoUsers), noRealization);
swco_WOA = zeros(length(NoUsers), noRealization);
po_WOA = zeros(length(NoUsers), noRealization);
su_WOA = zeros(length(NoUsers), noRealization);

for iN = 1:length(NoUsers)
    
    % Number of users
    N = NoUsers(iN);
    rho = ones(1,N);        % preferences of different mobile users
    % local computing capabilities
    f_l = f_user(1, 1:N);
    % channel gain
    channelGain = zeros(N, noRealization);
    for iReal = 1:noRealization
        [hA, dArray, eNB_position, UE_position] = channelModel(N, cellRadiusMax, cellRadiusMin, logNormalMean, logNormalDeviation);
        channelGain(:,iReal) = hA;
    end

    % The simulation runs for noRealization number of realizations
    for iReal = 1:noRealization
        % print the realization number
        % fprintf('Star the realization no. %d.\n', iReal);

        % create channel models
        gArray = channelGain(:,iReal);

        % ======== heuristic offloading decision algorithm (HODA) =========
        [ P_HODA, F_HODA, Z_l, Z_r_HODA, N_offloading, N_local ] = HODA( gArray, p_max, p_min, n0, f_l, f0, alpha, beta, S, lambda, kappa, W );
        po_HODA(iN, iReal) = length(N_offloading)/N;
        swco_HODA_l = sum(Z_l(N_local));
        swco_HODA_r = sum(Z_r_HODA(N_offloading));
        swco_HODA(iN, iReal) = swco_HODA_l + swco_HODA_r;
        % system utility
        Leader_pos_HODA = zeros(1,N); Leader_pos_HODA(N_offloading) = 1;
        T_l = beta./f_l;
        E_l = kappa.*beta.*(f_l).^2;
        Z_r = Inf(1,N);
        T_t = zeros(1,N);
        T_e = zeros(1,N);
        E_r = zeros(1,N);
        for j = 1:length(N_offloading)    
            T_t(N_offloading(j)) = alpha/(W*log2(1 + P_HODA(N_offloading(j))*gArray(N_offloading(j))/n0));   % Eq. (6)
            T_e(N_offloading(j)) = beta/F_HODA(N_offloading(j)); 
            E_r(N_offloading(j)) = P_HODA(N_offloading(j))*T_t(N_offloading(j))/zeta;
            Z_r(N_offloading(j)) = lambda(1)*(T_t(N_offloading(j)) + T_e(N_offloading(j))) + lambda(2)*E_r(N_offloading(j));
        end
        % the total completion time and overall energy consumption
        T_i = Leader_pos_HODA.*(T_e + T_t) + (1-Leader_pos_HODA).*T_l;
        E_i = Leader_pos_HODA.*E_r + (1-Leader_pos_HODA).*E_l;
        su_HODA(iN,iReal) = sum(lambda(1)*(T_l - T_i)./T_l + lambda(2)*(E_l - E_i)./E_l);

        % ================= Whale Optimization Algorithm ==================
        % Compute computation overhead due to local computing
        T_l = beta./f_l;
        E_l = kappa.*beta.*(f_l).^2;
        Z_l = lambda_t*T_l + lambda_e*E_l; 

        % Coefficients defined in Eq. (17)
        eta = rho.*lambda_t*alpha./(W.*T_l);
        gamma = rho.*lambda_e*alpha./(W.*E_l*zeta);
        tau = rho.*lambda_t;
        % criterion to stop the bi-section method
        varepsilon = 1e-6;  
        phi = @(y,a,x,eta) y*log2(1 + a*x) - (a/log(2))*(eta + y*x)/(1 + a*x);

        % Transmit power control 
        Nm = 1:N;
        P = zeros(1, N);
        for i = 1:length(Nm)
            if phi(gamma(Nm(i)), gArray(i), p_max, eta(Nm(i))) <= 0
                % line 3
                P(Nm(i)) = p_max;
            else
                % line 4
                p_s = p_min;    
                p_t = p_max;
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

        % % save the parameters and then call them in the Get_Functions_details.m 
        % save('BWOA_HODA_Paras','alpha','beta','f_l','kappa','lambda','W','gArray','P','zeta','S');

        [dim, fobj] = CostFunction_Ex3(N, S, gArray, P, n0, W, alpha, beta, f_l, kappa, lambda, zeta);
        [Leader_score,Leader_pos,Leader_pos_F,Convergence_curve] = BWOA_HODA(SearchAgents_no,Max_iter,BWOA_num,fobj,dim,rho,lambda,f_l,f0);
        Z_r = Inf(1,N);
        T_t = zeros(1,N);
        T_e = zeros(1,N);
        E_r = zeros(1,N);
        N_H = find(Leader_pos == 1);
        for j = 1:length(N_H)    
            T_t(N_H(j)) = alpha/(W*log2(1 + P(N_H(j))*gArray(N_H(j))/n0));   % Eq. (6)
            T_e(N_H(j)) = beta/Leader_pos_F(N_H(j)); 
            E_r(N_H(j)) = P(N_H(j))*T_t(N_H(j))/zeta;
            Z_r(N_H(j)) = lambda(1)*(T_t(N_H(j)) + T_e(N_H(j))) + lambda(2)*E_r(N_H(j));
        end
        % the total completion time and overall energy consumption
        T_i = Leader_pos.*(T_e + T_t) + (1-Leader_pos).*T_l;
        E_i = Leader_pos.*E_r + (1-Leader_pos).*E_l;
        % compute the outputs
        po_WOA(iN, iReal) = sum(Leader_pos > 0)/length(Leader_pos);
        swco_WOA_l = sum(Z_l(Leader_pos == 0));
        swco_WOA_r = sum(Z_r(Leader_pos == 1));
        swco_WOA(iN, iReal) = swco_WOA_l + swco_WOA_r;
        su_WOA(iN,iReal) = sum(lambda(1)*(T_l - T_i)./T_l + lambda(2)*(E_l - E_i)./E_l);
    end

end

% Finalize the performance metrics
po_HODA = mean(po_HODA,2);
su_HODA = mean(su_HODA,2);
swco_HODA = mean(swco_HODA,2);
po_WOA = mean(po_WOA,2);
su_WOA = mean(su_WOA,2);
swco_WOA = mean(swco_WOA,2);


% Plot the convergence evolution
% Plot the figure - the curve of the percentage of offloading users
i = 4;
figure(5*i + 1)
hold on
plot(1:length(po_WOA), po_WOA(1:length(po_WOA)), 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(1:length(po_HODA), po_HODA(1:length(po_HODA)), 'b-d', 'linewidth', 4.0, 'markers', 16);
set(gca,'FontSize',30,'XLim',[1 length(po_HODA)]);
xticks = 1:1:length(po_HODA);
set(gca,'xtick',xticks); 
% xticklabels({'2','4','6','8','10','12'})
xticklabels({'2','3','4','5','6','7','8','9','10'})
% set(gca,'xticklabel',xticklabels);
xlabel('Number of Users'); 
ylabel('Offloading Percentage');
legend('WOA-based Alg','HODA');
box on;
% Plot the figure - the curve of the system-wide computation overhead
figure(5*i + 2)
hold on
plot(1:length(swco_WOA), swco_WOA(1:length(swco_WOA)), 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(1:length(swco_HODA), swco_HODA(1:length(swco_HODA)), 'b-d', 'linewidth', 4.0, 'markers', 16);
set(gca,'FontSize',30,'XLim',[1 length(swco_HODA)]);
xticks = 1:1:length(swco_HODA);
set(gca,'xtick',xticks); 
xticklabels({'2','3','4','5','6','7','8','9','10'})
% set(gca,'xticklabel',xticklabels);
xlabel('Number of Users'); 
ylabel('Total Computation Overhead');
legend('WOA-based Alg','HODA');
box on;
% Plot the figure - the curve of the system utility
figure(5*i + 3)
hold on
plot(1:length(su_WOA), su_WOA(1:length(su_WOA)), 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(1:length(su_HODA), su_HODA(1:length(su_HODA)), 'b-d', 'linewidth', 4.0, 'markers', 16);
set(gca,'FontSize',30,'XLim',[1 length(su_HODA)]);
xticks = 1:1:length(su_HODA);
set(gca,'xtick',xticks); 
% xticklabels({'2','4','6','8','10','12'})
xticklabels({'2','3','4','5','6','7','8','9','10'})
set(gca,'xticklabel',xticklabels);
xlabel('Number of Users'); 
ylabel('System Utility');
legend('WOA-based Alg','HODA');
box on;


