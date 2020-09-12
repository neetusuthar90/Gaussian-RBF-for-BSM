clc;
close all;
clear all;
%err1 = [0.0004,0.0007,0.0003]; %MQ
err2 = [0.0316,0.0297,0.0184]; %American
err3 = [0.0798, 0.0737, 0.0612];%European
plot(err2,err3)
% plot(c1,err1,'-bo','MarkerFaceColor','b','MarkeredgeColor','b')
% hold on
% plot(c2,err2,'-bd','MarkerFaceColor','b','MarkeredgeColor','b')
% hold off
% box on
% ax = gca;
% ax.XGrid = 'off';
% ax.YGrid = 'on';
% xlabel('Shape parameter')
% ylabel('RMSE')
% legend('MQ-RBF','GA-RBF')