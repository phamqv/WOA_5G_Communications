%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %
%                                                                         %
%                                                                         %
%   Edited by: Quoc-Viet Pham                                             %
%   E-mail: vietpq90@gmail.com / vietpq@pusan.ac.kr                       %
%   Created Date: 4/1/2019                                                %
%   Current Date: 4/5/2019                                                %
%_________________________________________________________________________%


% The Whale Optimization Algorithm
function [Leader_score, Leader_pos, Convergence_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj)

% initialize position vector and score for the leader
Leader_pos = zeros(1, dim);
% change this to -inf for maximization problems
Leader_score = inf; 
Leader_pre_score = Leader_score;
% tolerance to stop the algorithm
delta = 1e-6;
Flag = 0;
% Initialize the positions of search agents. Size is SearchAgents_no x dim
Positions = initialization(SearchAgents_no, dim, ub, lb);
Convergence_curve = zeros(1,Max_iter);
% Loop counter
iter = 0;   

% Main loop
while iter < Max_iter && Flag <= 3
    for i = 1:size(Positions,1)

        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub = Positions(i,:) > ub;
        Flag4lb = Positions(i,:) < lb;
        Positions(i,:) = (Positions(i,:).*(~(Flag4ub+Flag4lb))) + ub.*Flag4ub + lb.*Flag4lb;
        
        % Calculate objective function for each search agent
        fitness = fobj(Positions(i,:));
        % fitness = fobj(Positions(i,:), I, gArray, Pc, varrho, n0, r_req);
        
        % Update the leader
        if fitness < Leader_score   % Change this to > for maximization problem
            Leader_score = fitness;         % Update alpha
            Leader_pos = Positions(i,:);
        end       
    end
    
    a = 2 - iter*((2)/Max_iter);    % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2 = -1 + iter*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i = 1:size(Positions,1)
        r1 = rand();    % r1 is a random number in [0,1]
        r2 = rand();    % r2 is a random number in [0,1]
        
        A = 2*a*r1-a;   % Eq. (2.3) in the paper
        C = 2*r2;       % Eq. (2.4) in the paper
        
        % parameters for spiral updating position
        b = 1;               %  parameters in Eq. (2.5)
        l = (a2-1)*rand + 1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j = 1:size(Positions,2)
            % follow the shrinking encircling mechanism or prey search
            if p < 0.5   
                % search for prey (exploration phase)
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j) - Positions(i,j)); % Eq. (2.7)
                    Positions(i,j) = X_rand(j) - A*D_X_rand;      % Eq. (2.8)
                % Shrinking encircling mechanism (exploitation phase)   
                elseif abs(A) < 1
                    D_Leader = abs(C*Leader_pos(j) - Positions(i,j)); % Eq. (2.1)
                    Positions(i,j) = Leader_pos(j) - A*D_Leader;      % Eq. (2.2)
                end
            % follow the spiral-shaped path
            elseif p>=0.5
                % Eq. (2.5)
                distance2Leader = abs(Leader_pos(j)-Positions(i,j));
                Positions(i,j) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + Leader_pos(j);
            end
        end
    end
    % increase the iteration index by 1
    iter = iter + 1;
    % negate the objective value (minimization --> maximization problem)
    Convergence_curve(iter) = -Leader_score;
    [iter -Leader_score]
    
    % check to see whether the stopping criterion is satisifed
    if abs(Leader_score - Leader_pre_score) < delta
        Flag = Flag + 1;
        Convergence_curve = Convergence_curve(1,1:iter);
        Leader_score = -Leader_score;
    else
        Leader_pre_score = Leader_score;
    end
end
% Leader_score = -Leader_score;

