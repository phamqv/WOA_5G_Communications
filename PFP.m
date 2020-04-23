function [p, Phi, iter] = PFP( I, pMax, pMin, n0, gArray, gArray_EV, stopEps)
% General Algorithm Procedure (GAP) for EE-SE tradeoff optimization problems
%
% =========================================================================================
% Related Journal Reference:
% 
% Zhichao Sheng et al., "Power Allocation for Energy Efficiency and Secrecy of Wireless Interference Networks," IEEE Transactions on Wireless Communications, vol. 17, no. 6, June 2018.
% Q.-V. Pham, W.-J. Hwang, "Fairness-Aware Spectral and Energy Efficiency in Spectrum-Sharing Wireless Networks," IEEE Transactions on Vehicular Technology, vol. 66, no. 11, pp. 10207-10219, Nov. 2017.
%
% All rights belong to Quoc-Viet Pham (email: vietpq90@gmail.com.com).
%
% This simulation code can be freely modified and distributed, as long as the the copyright notice
% section is kept unchanged. Any comments/suggestions are welcome. However, the original authors of 
% the code are not responsible for any damages caused by this code.
%
% If this simulation code (or its variation) is used in adademic research, one
% or more of the above journal/conference references should be appropriately cited.
%
% Author: Quoc-Viet Pham
% Personal site: https://sites.google.com/view/vietpq90
% Affiliation: Research Institute of Computer, Information and Communication 
%              Pusan National University
% Email: vietpq90@gmail.com / vietpq@pusan.ac.kr
% Created date: April 04, 2015
% Current date: April 08, 2019
% ========================================================================================

% ========================================================================================
% I: number of users
% W: baseband bandwidth
% pMin: maximum transmission power of each transmitter (assume same for every one)
% pMax: maximum transmission power of each transmitter (assume same for every one)
% distanceArray: the distance matrix between each transmitter and receiver (IxI)
% gridUnit: the length of a smallest grid
% attenFactor: the attenuation between two nodes is distance^(-attenFactor)
% hArray: path attenuation matrix (IxI) 
% h_ji is the channel gain from the transmitter j to the receiver i
% ========================================================================================

% Diagonal of channel gains matrix
diagGArray = reshape(diag(gArray),1,I);             

% energy efficiency
p = zeros(10000, I);
Phi = zeros(10000, 1);
Phi_pre = inf;

% initialization of the PFP algorithm
iter = 0;
Flag = 0;
% an initial feasible point
p_k = rand()*pMax;  

while Flag <= 3
    % Increase the iteration index by 1
    iter = iter + 1;
    
    % Solve the convex optimization problem (11) to obtain p^(k)
    x_k = diagGArray.*p_k./(p_k*gArray' - p_k.*diag(gArray)' + n0);
    xe_k = gArray_EV.*p_k./(p_k*gArray_EV' - p_k.*gArray_EV' + n0);
    % Eq. (8)
    f_k = @(p) log(1 + x_k) + (x_k./(1+ x_k)).*(2 - p_k./p - (p*gArray' - p.*diag(gArray)' + n0)./(p_k*gArray' - p_k.*diag(gArray)' + n0));
    % Eq. (10)
    g_k = @(p) log(1 + xe_k) + (1./(1+ xe_k)).*(0.5*gArray_EV.*((p.^2)./p_k + p_k)./(p*gArray_EV' - p.*gArray_EV + n0) - xe_k);
    % The function to be minimized
    fun = @(p) -(f_k(p) - g_k(p));
    lb = pMin;  % lower bound
    ub = pMax;  % upper bound
    [p_out] = fminimax(fun, p_k, [], [], [], [], lb, ub);
    
    % calculate the R_min^k as the value of the objective in (4) at p^(k)
    p(iter,:) = p_out;
    p_k = p_out;
%     Phi(iter) = -min(fval);
%     Phi_cur = Phi(iter);
    
    % data rate
    SINR = p_out.*diagGArray./(n0 + p_out*gArray' - p_out.*diagGArray);
    R = log(1 + SINR);
    %  wiredtapped data rate
    SINR_EV = p_out.*gArray_EV./(n0 + p_out*gArray_EV' - p_out.*gArray_EV);
    Gamma = log(1 + SINR_EV);
    % min secrecy throughput
    Phi(iter) = min(R - Gamma);
    Phi_cur = Phi(iter);
    
    
    % check to see whether the stopping criterion is satisfied
    if abs(Phi_cur - Phi_pre)/Phi_pre < stopEps
        Flag = Flag + 1;
    end
    Phi_pre = Phi_cur;
end    

% Return the output 
Phi = Phi(1:iter,:);
p = p(1:iter,:);

