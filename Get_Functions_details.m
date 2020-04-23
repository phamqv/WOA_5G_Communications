%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA)                                     %
%                                                                         %
%  Developed in MATLAB R2017a                                             %
%                                                                         %
%  Author and programmer: Quoc-Viet Pham                                  %
%                                                                         %
%         e-Mail: vietpq90@gmail.com                                      %
%                 vietpq@pusan.ac.kr                                      %
%                                                                         %
%       Homepage: https://sites.google.com/view/vietpq90/                 %
%                                                                         %
%   Main paper: Quoc-Viet Pham et al.                                     %
%               Whale Optimization Algorithm with Applications to         %
%               Resource Allocation in Wireless Networks,                 %
%               IEEE Transactions on Vehicular Technology, 2020           %
%               URL: https://ieeexplore.ieee.org/document/8993843         %
%                                                                         %
%_________________________________________________________________________%
%
% You are free to use the code provided you cite our IEEE TVT paper
% Q.-V. Pham, S. Mirjalili, N. Kumar, M. Alazab, and W.-J. Hwang, "Whale Optimization Algorithm with Applications to Resource Allocation in Wireless Networks," IEEE Transactions on Vehicular Technology, vol. 69, no. 4, pp. 4285-4297, Apr. 2020.
% bibtex: @ARTICLE{8993843, author={Q. {Pham} and S. {Mirjalili} and N. {Kumar} and M. {Alazab} and W. {Hwang}}, journal={IEEE Transactions on Vehicular Technology}, title={Whale Optimization Algorithm With Applications to Resource Allocation in Wireless Networks}, year={2020}, month={Apr.}, volume={69}, number={4}, pages={4285-4297},} 
%
%_________________________________________________________________________%

% This function containts full information and implementations of the benchmark 
% functions in Table 1, Table 2, and Table 3 in the paper

% lb is the lower bound: lb=[lb_1,lb_2,...,lb_d]
% up is the uppper bound: ub=[ub_1,ub_2,...,ub_d]
% dim is the number of variables (dimension of the problem)
% gArray_EV: [1 x N], where N is the number of users 

function [lb, ub, dim, fobj] = Get_Functions_details(FN, I, pMax, pMin, gArray, Pc, varrho, n0, r_req, gArray_EV)
 
    switch FN
        case 'Ex1'
            fobj = @F_Ex1;
            lb = pMin;
            ub = pMax;
            dim = I;
            
        case 'Ex2'
            % return the output
            fobj = @F_Ex2;
            lb = pMin;
            ub = pMax;
            dim = I;
    end
    
    % objective function - example 1
    function objf = F_Ex1(P)
        diagGArray = reshape(diag(gArray),1,I);
        % data rate
        SINR = P.*diagGArray./(n0 + P*gArray' - P.*diagGArray);
        R = log(1 + SINR);
        %  wiredtapped data rate
        SINR_EV = P.*gArray_EV./(n0 + P*gArray_EV' - P.*gArray_EV);
        Gamma = log(1 + SINR_EV);
        % secrecy data rate
        Phi = R - Gamma;

        % objective function (maximization --> minimization problem)
        % unconstrained problem --> no penalty factor
        objf = max(-Phi);
    end
    
    % objective function - example 2
    function objf = F_Ex2(P)
        diagGArray = reshape(diag(gArray),1,I);
        gamma = P.*diagGArray./(n0 + P*gArray' - P.*diagGArray);

        % penalty method to deal with inequality constraints
        % X.-S. Yang Nature-Inspired Optimization Algorithms (2014, Elsevier)
        muy = 1e14;     % muy can be taken as 10^13 to 10^15 
        H = (r_req > log2(1 + gamma));
        penalty = muy*sum(H.*((r_req - log2(1 + gamma)).^2));

        % objective function (maximization --> minimization problem)
        objf = -sum(log2(1 + gamma))/sum(Pc + varrho*P) + penalty;
    end
end