%_________________________________________________________________________%
%  Binary Whale Optimization Algorithm (WOA) source codes demo 1.0        %
%                                                                         %
%  Developed in MATLAB R2017a                                             %
%                                                                         %
%  Author and programmer: Quoc-Viet Pham                                  %
%                                                                         %
%         e-Mail: vietpq90@gmail.com                                      %
%                 vietpq90@changwon.ac.kr                                 %
%                                                                         %
%       Homepage: https://sites.google.com/view/vietpq90/                 %
%                                                                         %
%   Main paper: Quoc-Viet Pham                                            %
%               Whale Optimization Algorithm with Applications to         %
%               Resource Allocation in Wireless Networks,                 %
%               To be submitted to an IEEE Journal                        %
%               DOI:                                                      %
%                                                                         %
%_________________________________________________________________________%


% The Whale Optimization Algorithm
function [Leader_score,Leader_pos,Leader_pos_F,Convergence_curve] = BWOA_HODA(SearchAgents_no,Max_iter,BWOA_num,fobj,dim,rho,lambda,f_l,f0)

% initialize position vector and score for the leader
Leader_pos = zeros(1,dim);
Leader_pos_F = zeros(1,dim);
% change this to -inf for maximization problems
Leader_score = inf; 
Leader_score_pre = Leader_score;
% tolerance to stop the algorithm
delta = 1e-6;
todoTol = 0;
Flag = 0;
% Initialize the positions of search agents.
Positions = zeros(SearchAgents_no,dim);
for i = 1:size(Positions,1) 	% For each seach agent
    for j = 1:size(Positions,2) % For each variable
        if rand <= 0.5
            Positions(i,j) = 0;
        else
            Positions(i,j) = 1;
        end
    end
end
Convergence_curve = zeros(1,Max_iter);
C_Step = zeros(SearchAgents_no,dim);
% Loop counter
iter = 0;
tau = rho.*lambda(1);
% Main loop
while iter < Max_iter && Flag <= 3
	% Update the computing resource allocation
    Pos_F = zeros(SearchAgents_no,dim);
    for i = 1:size(Positions,1)
        for j = 1:size(Positions,2)
            if Positions(i,j) == 1
                Pos_F(i,j) = f0.*sqrt(tau(j).*f_l(j))./sum(sqrt(tau.*f_l.*Positions(i,:)));
            end
        end
    end

    Tmp = zeros(1,size(Positions,1));
    for i = 1:size(Positions,1)
        % Calculate objective function for each search agent
        Pos = [Positions(i,:),Pos_F(i,:)];
        fitness = fobj(Pos);
        Tmp(i) = fitness;
        
        % Update the leader
        if fitness < Leader_score % Change this to > for maximization problem
            Leader_score = fitness;         % Update alpha
            Leader_pos = Positions(i,:);
            Leader_pos_F = Pos_F(i,:);
        end       
    end
    
    a = 2-iter*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2 = -1+iter*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i = 1:size(Positions,1)
        r1 = rand();    % r1 is a random number in [0,1]
        r2 = rand();    % r2 is a random number in [0,1]
        
        A = 2*a*r1-a;   % Eq. (2.3) in the paper
        C = 2*r2;       % Eq. (2.4) in the paper
        
        % parameters for spiral updating position
        b = 1;               	%  parameters in Eq. (2.5)
        l = (a2-1)*rand + 1;   	%  parameters in Eq. (2.5)
        
        p = rand();        		% p in Eq. (2.6)
        
        for j = 1:size(Positions,2)
            % follow the shrinking encircling mechanism or prey search
            if p < 0.5   
                % search for prey (exploration phase)
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j) - Positions(i,j));
					% C_Step(i,j) = D_X_rand;
                    C_Step(i,j) = X_rand(j) - A*D_X_rand;
                % Shrinking encircling mechanism (exploitation phase)   
                elseif abs(A) < 1
                    D_Leader = abs(C*Leader_pos(j) - Positions(i,j)); 	% Eq. (2.1)
					% C_Step(i,j) = D_Leader;
                    C_Step(i,j) = Leader_pos(j) - A*D_Leader;
                end
            % follow the spiral-shaped path (exploitation phase)
            elseif p >= 0.5
                distance2Leader = abs(Leader_pos(j)-Positions(i,j));	% Eq. (2.5)	
				% C_Step(i,j) = distance2Leader;
                C_Step(i,j) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + Leader_pos(j);
            end
			
            if BWOA_num == 1
                % s = 1/(1+exp(-2*A*C_Step(i,j)));			% S1 transfer function           
                s = 1/(1+exp(-2*C_Step(i,j)));
            end
            if BWOA_num == 2
                % s = 1/(1+exp(-A*C_Step(i,j)));   			% S2 transfer function              
                s = 1/(1+exp(-C_Step(i,j)));
            end
            if BWOA_num == 3
                % s = 1/(1+exp(-A*C_Step(i,j)/2)); 			% S3 transfer function              
                s = 1/(1+exp(-C_Step(i,j)/2));
            end
            if BWOA_num == 4
                % s = 1/(1+exp(-A*C_Step(i,j)/3));  			% S4 transfer function
                s = 1/(1+exp(-C_Step(i,j)/3));
            end			 			
            % S-shaped transfer functions
            if BWOA_num <= 4 
                if rand < s 	% Check Subsection II-D. Binary Whale Optimization Algorithm
                    Positions(i,j) = 1;
                else
                    Positions(i,j) = 0;
                end
            end
			 
            if BWOA_num == 5
                % s = abs(erf(((sqrt(pi)/2)*A*C_Step(i,j)))); 	% V1 transfer function
                s = abs(erf(((sqrt(pi)/2)*C_Step(i,j))));
            end 
            if BWOA_num == 6
                % s = abs(tanh(A*C_Step(i,j))); 					% V2 transfer function
                s = abs(tanh(C_Step(i,j)));
            end            
            if BWOA_num == 7
                % s = abs(A*C_Step(i,j)/sqrt((1+A*C_Step(i,j)^2))); % V3 transfer function
                s = abs(C_Step(i,j)/sqrt((1+C_Step(i,j)^2)));
            end            
            if BWOA_num == 8
                % s = abs((2/pi)*atan((pi/2)*A*C_Step(i,j)));    % V4 transfer function (VPSO)         
                s = abs((2/pi)*atan((pi/2)*C_Step(i,j)));
            end
            if BWOA_num == 9
                % s = 1/(1+exp(-10*(A*C_Step(i,j)-0.5)));	     
                s = 1/(1+exp(-10*(C_Step(i,j)-0.5)));
            end 
            % V-shaped transfer functions
            p_rand = rand();
            if BWOA_num > 4 && BWOA_num <= 8 || BWOA_num == 9	
                if p_rand < s 
                    Positions(i,j) = ~Positions(i,j); 
                end
            end
        end
    end
    % increase the iteration index by 1
    iter = iter + 1;
    Convergence_curve(iter) = Leader_score;
    [iter Leader_score]
    
    if todoTol == 1 && abs(Leader_score - Leader_score_pre) < delta
        Flag = Flag + 1;
		Convergence_curve = Convergence_curve(1,1:iter);
    end
    Leader_score_pre = Leader_score;
end

end



