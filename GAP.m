function [P, eta, fhp, gamma, iterNum_inner, iterNum_outer] = GAP( I, pMin, pMax, Pc, r_req, n0, varrho, gArray, stopEps)
% General Algorithm Procedure (GAP) for EE-SE tradeoff optimization problems
%
% =========================================================================================
% Related Journal Reference:
% 
% Y. Li, M. Sheng, C. Yang, and X. Wang, "Energy efficiency and spectral efficiency tradeoff in interference-limited wireless networks," Communications Letters, IEEE, vol. 17, no. 10, pp. 1924–1927, Oct 2013.
% Q.-V. Pham, W.-J. Hwang, "Fairness-Aware Spectral and Energy Efficiency in Spectrum-Sharing Wireless Networks," IEEE Transactions on Vehicular Technology, vol. 66, no. 11, pp. 10207-10219, Nov. 2017.
%
% All rights belong to Quoc-Viet Pham (email: vietpq90@gmail.com).
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
% h_ji is the channel gain from the transmitter i to the receiver j
% ========================================================================================

% Diagonal of channel gains matrix
diagGArray = reshape(diag(gArray),1,I);             

% energy efficiency
eta = zeros(10000, 1);
% the maximum tolerance
delta = stopEps;        

% initializatin of the GAP algorithm
iter_Tighten = 0;
FLAG = 0;
% maximum and minimum energy efficiency, eta_max and eta_min, respectively
% if the maximum power is defined as a constaint
etaMax = sum(log2(1 + (pMax./(n0 + pMin*gArray' - pMin.*diagGArray)).*diagGArray))/sum(Pc + varrho*pMin.*ones(1, I));
etaMin = sum(log2(1 + (pMin./(n0 + pMax*gArray' - pMax.*diagGArray)).*diagGArray))/sum(Pc + varrho*pMax.*ones(1, I));

% defined functions of the D.C. approximation method
fP = @(P, Pc, G, n0, eta, varrho) log2(exp(1))*sum(log(n0 + P*G')) - eta*sum(varrho*P + Pc);
hP = @(P, G, n0) log2(exp(1))*sum(log(n0 + P*G' - P.*diag(G)'));
% e_i is an M-dimensional column vector
e = zeros(I, I);
for i = 1:I
    e(i,:) = gArray(i,:)/log(2);
    e(i,i) = 0;
end

while (FLAG < 1)
    % iteration
    iter_Tighten = iter_Tighten + 1;
    iter = 0;
    eta(iter_Tighten, :) = (etaMax + etaMin)/2;
    eta_cur = eta(iter_Tighten, :);
    
    % Initialization
    P = zeros(10000, I);        % power and outage prices
    P_cur = pMax.*rand(1, I);   % current transmit power
    gamma = zeros(10000, I);    % SINR
    gamma_cur = zeros(1, I);
    fhp = zeros(1000, 1);       % substraction of (fp-hp)
    I_gap_cur = fP(P_cur, Pc, gArray, n0, eta_cur, varrho) - hP(P_cur, gArray, n0);
    I_gap_pre = I_gap_cur;
    FLAG_IN = 0;   
    while (FLAG_IN == 0)
        % iteration
        iter = iter + 1;
        
        % convex software packages
        % step 3 - solve problem (16) to obtain the solution P*
        cvx_begin
            variable P_cvx(1, I)
            maximize (log2(exp(1))*sum(log(n0 + P_cvx*gArray')) - eta_cur*sum(varrho*P_cvx + Pc)...
                - (log2(exp(1))*sum(log(n0 + P_cur*gArray' - P_cur.*diagGArray)) + ((1./(n0 + P_cur*gArray' - P_cur.*diagGArray))*e)*(P_cvx - P_cur)'))
            subject to 
                P_cvx.*diagGArray + (1 - 2.^r_req).*(n0 + P_cvx*gArray' - P_cvx.*diagGArray) >= 0;
                P_cvx >= pMin; 
                P_cvx <= pMax;
        cvx_end
        
        % assign the optimal solution of cvx solver to the current power
        P(iter, :) = P_cvx;
        P_cur = P(iter,:);
        
        % compute SINR
        gamma(iter,:) = P_cur.*diagGArray./(n0 + P_cur*gArray' - P_cur.*diagGArray);
        gamma_cur = gamma(iter,:);
        
        % update and check the stopping condition
        I_gap_cur = abs(fP(P_cur, Pc, gArray, n0, eta_cur, varrho) - hP(P_cur, gArray, n0));
        fhp(iter,:) = fP(P_cur, Pc, gArray, n0, eta_cur, varrho) - hP(P_cur, gArray, n0);
        if (abs(I_gap_cur - I_gap_pre) < delta)
            FLAG_IN = 1;
        else
            I_gap_pre = I_gap_cur;
        end
    end    % end of the Iterative Power Allocation Algorithm (IPAA)
    iterNum_inner = iter;
    
    m_gap = sum(log2(1 + gamma_cur)) - eta_cur*sum(Pc + varrho*P_cur);
    % update the energy efficiency
    if (abs(m_gap) <= stopEps)
        P = P(1:iterNum_inner,:);
        gamma = gamma(1:iterNum_inner,:);
        fhp = fhp(1:iterNum_inner,:);
        eta(iter_Tighten, :) = sum(log2(1 + gamma_cur))/sum(Pc + varrho*P_cur);
        FLAG = 2;
        iter_Tighten_Out = iter_Tighten;
    else
        if (m_gap < 0)
            etaMax = eta_cur;
        else
            etaMin = eta_cur;
        end
    end
end    % end of the General Algorithm Procedure (GAP) for EE-SE tradeoff optimization problems
% number of outer loop
iterNum_outer = iter_Tighten_Out;
eta = eta(1:iter_Tighten_Out,:);

