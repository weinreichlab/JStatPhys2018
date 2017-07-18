
% This code plots the observed values tau_b and p-value for the
% Palmer dataset as analyzed by me, and the corresponding Tau_b c.d.f. 
%
% This is how I produced the plot in Fig 1c.
%
% DMW June 20, 2017
%
% Make plot lines wider (2 point) and increase text point size (24 point)
% in order to make legible after shrinking
%
% June 23, 2017

close all
clear all

% These first two values comes from TauTable.csv, which is produced by
% BuildKendallTauTable.m
ObsTaub = 0.1921;
SizePalmerRanks = 55;
genotypes = 64;
loci = log2(genotypes);

% Now we need the null distribution of Tau_b values, which may already
% exist on the disk... Here's the file name convention used in
% BuildKendallTauTable.m
NullTaubFilename = ...
    sprintf('./NullTau_bDistrs/%dLoci%dSigWalsh.mat',loci,...
    SizePalmerRanks);
% First see if directory exists, and create it if necessary
if ~exist('./NullTau_bDistrs')
    mkdir('./NullTau_bDistrs');
end
% Now see if this particular null distribution has already been
% simulated
if exist(NullTaubFilename)
    % if so, use it
    load(NullTaubFilename)
    total_pseudo_replicates = size(taub_null_distr,1);
else
    % Build the canonical unitation for this dataset
    canonical_unitations = [];
    for i=0:loci
        for j=1:nchoosek(loci,i)
            canonical_unitations(end+1) = i;
        end
    end
    % Now build the null distribution
    total_pseudo_replicates = 1e5;
    taub_null_distr = zeros(total_pseudo_replicates,1);
    for i = 1:total_pseudo_replicates
        test_sample = randsample(canonical_unitations, SizePalmerRanks);
        test_sample2 = randsample(canonical_unitations, SizePalmerRanks);
        taub_null_distr(i) = ktaub([test_sample' test_sample2'],0.05,0); 
    end
    taub_null_distr = sort(taub_null_distr);
    % And save it for next time
    save(NullTaubFilename,'taub_null_distr');
end

% Finally we can get P
for i= 1:total_pseudo_replicates
    if ObsTaub < taub_null_distr(i)
        break
    end
end
% not really: at this point pval is just integer index to plot with, 
% rather than the (0,1) probability itself.
pval1 = total_pseudo_replicates - i;

% plot([1:-1/total_pseudo_replicates:1/total_pseudo_replicates],...
plot([total_pseudo_replicates:-1:1],taub_null_distr,'--','LineWidth',2)
hold on
plot([pval1,pval1],[-.5,ObsTaub],'r','LineWidth',2)
plot([0,pval1],[ObsTaub,ObsTaub],'r','LineWidth',2)
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