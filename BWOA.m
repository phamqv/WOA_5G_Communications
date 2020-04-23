%_________________________________________________________________________%
%  Binary Whale Optimization Algorithm (WOA)                              %
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

% The Whale Optimization Algorithm
function [Leader_score,Leader_pos,Convergence_curve] = BWOA(SearchAgents_no,Max_iter,BWOA_num,fobj,dim)

% initialize position vector and score for the leader
Leader_pos = zeros(1,dim);
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
D_X = zeros(SearchAgents_no,dim);
% Loop counter
iter = 0;

% Main loop
while iter < Max_iter && Flag <= 3
    for i = 1:size(Positions,1)

        % Calculate objective function for each search agent
        fitness = fobj(Positions(i,:));
        
        % Update the leader
        if fitness < Leader_score   % Change this to > for maximization problem
            Leader_score = fitness;         % Update the best score
            Leader_pos = Positions(i,:);    % Update the position
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
        
        p = rand();        		% p in Eq. (2.6)
        
        for j = 1:size(Positions,2)
            % follow the shrinking encircling mechanism or prey search
            if p < 0.5   
                % search for prey (exploration phase)
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j) - Positions(i,j)); 		% Eq. (2.7)
					D_X(i,j) = D_X_rand;
                % Shrinking encircling mechanism (exploitation phase)   
                elseif abs(A) < 1
                    D_Leader = abs(C*Leader_pos(j) - Positions(i,j)); 	% Eq. (2.1)
					D_X(i,j) = D_Leader;
                end
            % follow the spiral-shaped path (exploration phase)
            elseif p >= 0.5
                distance2Leader = abs(Leader_pos(j)-Positions(i,j));	% Eq. (2.5)	
				D_X(i,j) = distance2Leader;
            end
			
            if BWOA_num == 1
                s = 1/(1+exp(-2*A*D_X(i,j)));			% S1 transfer function           
            end
            if BWOA_num == 2
                s = 1/(1+exp(-A*D_X(i,j)));   			% S2 transfer function              
            end
            if BWOA_num == 3
                s = 1/(1+exp(-A*D_X(i,j)/2)); 			% S3 transfer function              
            end
            if BWOA_num == 4
               s = 1/(1+exp(-A*D_X(i,j)/3));  			% S4 transfer function
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
                s = abs(erf(((sqrt(pi)/2)*A*D_X(i,j)))); 	% V1 transfer function
            end 
            if BWOA_num == 6
                s = abs(tanh(A*D_X(i,j))); 					% V2 transfer function
            end            
            if BWOA_num == 7
                s = abs(A*D_X(i,j)/sqrt((1+A*D_X(i,j)^2))); % V3 transfer function
            end            
            if BWOA_num == 8
                s = abs((2/pi)*atan((pi/2)*A*D_X(i,j)));    % V4 transfer function (VPSO)         
            end
            if BWOA_num == 9
                s = 1/(1+exp(-10*(A*D_X(i,j)-0.5)));	     
            end 
            % V-shaped transfer functions
            if BWOA_num > 4 && BWOA_num <= 9 
                if rand < s 
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



