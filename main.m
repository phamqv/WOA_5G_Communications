%_________________________________________________________________________%
%  Main file                                                              %
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

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb = [lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub = [ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run WOA: [Best_score,Best_pos,WOA_cg_curve] = WOA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________________________________________

clear all 
clc

% parameters for the WOA algorithm
SearchAgents_no = 30;   % Number of search agents
Max_iter = 500;         % Maximum number of iterations

% parameters for simulations
I = 4;                          % number of users
Pc = 0.1*1e-3;                  % circuit power consumption
varrho = 1;                     % power-amplifier inefficiency factor
n0 = 0.1*1e-6;                  % noise power
pMax = 1e-3*[0.7 0.8 0.9 1];    % max transmit power
pMin = zeros(1, I);             % min transmit power
stopEps = 1e-3;                 % stopping criterion
r_req = 0.8;                    % minimum rate requirement
gArray = [0.4310 0.0002 0.0129 0.0011;
		  0.0002 0.3018 0.0005 0.0031;
		  0.2605 0.0008 0.4266 0.0099;
		  0.0039 0.0054 0.1007 0.0634];
diagGArray = reshape(diag(gArray),1,I);
% d_EV = 10*rand(1,I). Random links are created in a 10m-by-10m area
d_EV = [4.1454, 3.6180, 6.4587, 4.9546];
gArray_EV = d_EV.^-4;

% =========================================================================
% =============================== EXAMPLE 1 ===============================
% =========================================================================

% Load details of the selected benchmark function for the WOA algorithm
Function_name = 'Ex1';
p_Max = 10e-3*ones(1,I);
p_Min = 1e-10*ones(1,I);
% WOA-based Algorithm
[lb, ub, dim, fobj] = Get_Functions_details(Function_name, I, p_Max, p_Min, gArray, Pc, varrho, n0, r_req, gArray_EV);
[Best_score, Best_pos, WOA_cg_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
display(['The best solution obtained by WOA is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by WOA is : ', num2str(Best_score)]);
% Path-following procedure (PFP)
[~, Phi, ~] = PFP( I, p_Max, p_Min, n0, gArray, gArray_EV, stopEps);
% ============================== plot figures =============================
i = 0;
% Plot the figure - the curve of min secrecy throughput
figure(5*i + 1)
hold on
plot(0:length(WOA_cg_curve), [0 WOA_cg_curve], 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(0:length(Phi), [0 Phi'], 'b-d', 'linewidth', 4.0, 'markers', 16);
xticks = 0:2:length(WOA_cg_curve);
set(gca,'xtick',xticks); 
set(gca,'FontSize',30,'XLim',[0 length(WOA_cg_curve)]);
xlabel('Iteration'); 
ylabel('Min secrecy throughput(bps/Hz)');
legend('WOA-based Alg','Path-Following Procedure');
box on;



% =========================================================================
% =============================== EXAMPLE 2 ===============================
% =========================================================================

% Load details of the selected benchmark function for the WOA algorithm
Function_name = 'Ex2';
[lb, ub, dim, fobj] = Get_Functions_details(Function_name, I, pMax, pMin, gArray, Pc, varrho, n0, r_req, gArray_EV);
[Best_score, Best_pos, WOA_cg_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);

display(['The best solution obtained by WOA is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by WOA is : ', num2str(Best_score)]);
% ============================== plot figures =============================
i = 1;
% Plot the figure - the curve of energy efficiency
figure(5*i + 1)
semilogy(WOA_cg_curve,'Color','r')
title('Objective space')
xlabel('Iteration');
ylabel('Best GEE obtained so far');
axis tight
grid on
box on
legend('WOA')

% The convergence evolution of the WOA-based algorithm
R_rq = [0.0 0.5 1.0];
WOA_curve = cell(length(R_rq),1);
Function_name = 'Ex2';
for i = 1:length(R_rq)
    r_req = R_rq(i);
    
    % Load the objective and call the WOA-based algorithm
    [lb, ub, dim, fobj] = Get_Functions_details(Function_name, I, pMax, pMin, gArray, Pc, varrho, n0, r_req, gArray_EV);
    [Best_score, Best_pos, WOA_cg_curve] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
    WOA_curve{i} = WOA_cg_curve;
end
% Plot the convergence evolution
% Plot the figure - the curve of energy efficiency
figure(5*i + 2)
hold on
plot(0:length(WOA_curve{1}), [0 WOA_curve{1}(1:length(WOA_curve{1}))/1e4], 'r-p', 'linewidth', 3.0, 'markers', 10);
plot(0:length(WOA_curve{2}), [0 WOA_curve{2}(1:length(WOA_curve{2}))/1e4], 'b-d', 'linewidth', 3.0, 'markers', 10);
plot(0:length(WOA_curve{3}), [0 WOA_curve{3}(1:length(WOA_curve{3}))/1e4], 'm-o', 'linewidth', 3.0, 'markers', 10);
% plot(1:length(WOA_curve{1}), WOA_curve{1}(1:length(WOA_curve{1})), 'r-s', 'linewidth', 4.0);
% plot(1:length(WOA_curve{2}), WOA_curve{2}(1:length(WOA_curve{2})), 'g-d', 'linewidth', 4.0);
% plot(1:length(WOA_curve{3}), WOA_curve{3}(1:length(WOA_curve{3})), 'm-o', 'linewidth', 4.0);
set(gca,'FontSize',30,'XLim',[0 48]);
xticks = 0:6:length(WOA_cg_curve);
set(gca,'xtick',xticks); 
xlabel('Iteration'); 
ylabel('Best GEE obtained so far (x10^4)');
legend('R_i^{req} = 0.0 bit/s/Hz','R_i^{req} = 0.5 bit/s/Hz','R_i^{req} = 1.0 bit/s/Hz');
box on;

% The convergence evolution of the WOA-based algorithm
R_rq = 0.0:0.2:1.0;
WOA_curve = zeros(length(R_rq),1);
GAP_curve = zeros(length(R_rq),1);
Function_name = 'Ex2';
for i = 1:length(R_rq)
    r_req = R_rq(i);
    
    % Load the objective and call the WOA-based algorithm
    [lb, ub, dim, fobj] = Get_Functions_details(Function_name, I, pMax, pMin, gArray, Pc, varrho, n0, r_req, gArray_EV);
    [Best_score, ~, ~] = WOA(SearchAgents_no, Max_iter, lb, ub, dim, fobj);
    WOA_curve(i) = Best_score;
    
    % Call the GAP algroithm from [2]
    [~, eta, ~, ~] = GAP( I, pMin, pMax, Pc, r_req, n0, varrho, gArray, stopEps);
    GAP_curve(i) = eta(end);
end
% Plot the convergence evolution
% Plot the figure - the curve of energy efficiency
figure(5*i + 3)
hold on
plot(1:length(WOA_curve), WOA_curve(1:length(WOA_curve))/1e4, 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(1:length(GAP_curve), GAP_curve(1:length(GAP_curve))/1e4, 'b-d', 'linewidth', 4.0, 'markers', 16);
set(gca,'FontSize',30,'XLim',[1 length(R_rq)]);
xticks = 1:1:length(R_rq);
set(gca,'xtick',xticks); 
xticklabels({'0','0.2','0.4','0.6','0.8','1.0'})
set(gca,'xticklabel',xticklabels);
xlabel('Minimum Rate (bit/s/Hz)'); 
ylabel('Energy Efficiency (bit/Joule/Hz)');
legend('WOA-based Alg','GAP');
box on;
