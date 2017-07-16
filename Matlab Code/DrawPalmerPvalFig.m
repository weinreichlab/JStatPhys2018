
% This code plots the observed values tau_b and p-value for the
% Palmer dataset and the corresponding Tau_b c.d.f. 
%
% DMW June 20, 2017
%
% Make plot lines wider (2 point) and increase text point size (24 point)
% in order to make legible after shrinking
%
% June 23, 2017

close all

PALMER_TAU_B = 0.1980;
PALMER_P = 0.03787;
PALMER_SIG_COEFFICIENTS = 54;

load('./NullTau_bDistr.mat');
total_pseudo_replicates = size(taub_null_distr,1);
PALMER_P = PALMER_P * total_pseudo_replicates;

% plot([1:-1/total_pseudo_replicates:1/total_pseudo_replicates],...
plot([total_pseudo_replicates:-1:1],...
    taub_null_distr(:,PALMER_SIG_COEFFICIENTS),'--','LineWidth',2)
hold on
plot([PALMER_P,PALMER_P],[-.5,PALMER_TAU_B],'r','LineWidth',2)
plot([0,PALMER_P],[PALMER_TAU_B,PALMER_TAU_B],'r','LineWidth',2)
set(gca, 'FontSize',24);
%xticks([1:total_pseudo_replicates/5:total_pseudo_replicates]);
yticks([-0.5:.25:.5]);

% Create xlabel and ylabel
X = xlabel({'Sorted Test Statistic Index'});
Y = ylabel({'Sorted Test Statistic Value'});
set(X, 'FontSize', 24, 'FontName','Helvetica');
set(Y, 'FontSize', 24, 'FontName','Helvetica');
%set(gca,'xscale','log')


% gcf is some kind of magic handle to the figure; I found it in the
% saveas help page w/o much explanation.
saveas (gcf,'../PalmerPvalFig','pdf');

hold off