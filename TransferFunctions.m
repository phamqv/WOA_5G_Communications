%_________________________________________________________________________%
%  Transfer functions                                                     %
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

x = -4:0.5:4;
T_1 = 1./(1 + exp(-2*x));
T_2 = 1./(1 + exp(-x));
T_3 = 1./(1 + exp(-x/2));

figure(1)
hold on
plot(1:length(x), T_1(1:length(x)), 'r-p', 'linewidth', 4.0, 'markers', 16);
plot(1:length(x), T_2(1:length(x)), 'b-d', 'linewidth', 4.0, 'markers', 16);
plot(1:length(x), T_3(1:length(x)), 'm-s', 'linewidth', 4.0, 'markers', 16);
set(gca,'FontSize',30,'XLim',[1 length(x)]);
xticks = 1:2:length(x);
set(gca,'xtick',xticks); 
xticklabels({'-4','-3','-2','-1','0','1','2','3','4'})
set(gca,'xticklabel',xticklabels);
xlabel('The distance x'); 
ylabel('The probability of changing T(x)');
legend('T2','T1','T1/2');
box on;