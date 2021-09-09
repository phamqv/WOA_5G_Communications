%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA)                                     %
%                                                                         %
%  Developed in MATLAB R2017a                                             %
%                                                                         %
%  Author and programmer: Quoc-Viet Pham                                  %
%                                                                         %
%       E-Mail: vieptq90@gmail.com                                        %
%                                                                         %
%       Homepage: https://sites.google.com/view/vietpq90                  %
%                                                                         %
%   Main paper: Quoc-Viet Pham and Won-Joo Hwang                          %
%               Title: WOA for EE-SE in Wireless Interference Networks    %
%               Journal:                                                  %
%               DOI:                                                      %
%                                                                         %
%_________________________________________________________________________%

function [dim, fobj] = CostFunction_Ex3(N, S, gArray, P, n0, W, alpha, beta, f_l, kappa, lambda, zeta)
 
dim = N;
fobj = @F_Ex3;

% objective function - example 3
function objf = F_Ex3(Pos)
    A = Pos(1:N);
    F = Pos(N+1:end);

    %===========Compute computation overhead due to local computing============
    T_l = beta./f_l;
    E_l = kappa.*beta.*(f_l).^2;
    Z_l = lambda(1)*T_l + lambda(2)*E_l;

    %===========Compute computation overhead due to local computing============
    Z_r = zeros(1,N);
    T_t = zeros(1,N);
    T_e = zeros(1,N);
    E_r = zeros(1,N);
    for j = 1:N
        if A(j) == 1
            % Eq. (6) - uplink transmission time
            T_t(j) = alpha/(W*log2(1 + P(j)*gArray(j)/n0)); 
            % Eq. (7) - remote completion time
            T_e(j) = beta/F(j); 	
            % Eq. (8) - Energy consumption for remote computation
            E_r(j) = (P(j)/zeta)*T_t(j);  	
            Z_r(j) = lambda(1)*(T_t(j) + T_e(j)) + lambda(2)*E_r(j);
        end
    end
    T_i = A.*(T_e + T_t) + (1-A).*T_l;
    E_i = A.*E_r + (1-A).*E_l;

    % penalty method to deal with inequality constraints
    % X.-S. Yang Nature-Inspired Optimization Algorithms (2014, Elsevier)
    muy = 1e14;     % muy can be taken as 10^13 to 10^15 
    H = (sum(A) > S);
    penalty = muy*sum(H.*((sum(A) - S).^2));

    % objective function (maximization --> minimization problem)
    % objf = A*Z_r' + (1-A)*Z_l' + penalty;
    objf = -sum(lambda(1)*(T_l - T_i)./T_l + lambda(2)*(E_l - E_i)./E_l) + penalty;
end

end

