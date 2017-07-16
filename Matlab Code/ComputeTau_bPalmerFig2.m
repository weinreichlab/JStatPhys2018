% This code is a modified form of our standard analysis pipeline built to
% allow analysis of the data in Fig 2a of Palmer et al. Those authors
% aggregated all epistatic effects of third-order and above. Here I've
% hard-coded their observations (PalmerRanks) and then modified the null
% expectation (canonical_unitations) to pool all ranks of 3 and above into
% 3's. Finally it does the requisite permutation test right here, since
% this null isn't any of the ones captured for the analysis of other data
% sets.
%
% DMW June 29, 2017

clear all;
close all;

canonical_unitations = [];
total_pseudo_replicates = 1e5; 
loci = 6;
genotypes = 96;

% These come directly from Fig 2a in Palmer
PalmerRanks = [0,1,2,3,3,3,3,3,1,2,3,2,3,1,3,3,2,2,1,3,3,3,2,...
    3,3,3,3,3,3,3,3,3,3,2,3,2,2,3,3,3,3,1,2,2,3,2,3,3,3,3,2,...
    3,3,3,2,3,3,3,3,3,3,1,3,3,3,2,1,3,3,3];
SizePalmerRanks = size(PalmerRanks,2);

% Build the canonical unitation for the Palmer-sized data. First the three
% unitations that are given in the figure...
for i=0:2
    for j=1:nchoosek(loci,i)
        canonical_unitations(end+1) = i;
    end
end
% Everything else is unitation 3.
i = 1;
for j=0:2
    i = i + nchoosek(loci,j);
end
for i=i:genotypes
    canonical_unitations(i) = 3;
end

% Compute observed taub
ObsTaub = ktaub([canonical_unitations(1:SizePalmerRanks)' ...
    PalmerRanks'],0.05,0);

% Now get the distribution
NullDistr = zeros(total_pseudo_replicates,1);
for i = 1:total_pseudo_replicates
    test_sample = randsample(canonical_unitations, SizePalmerRanks);
    test_sample2 = randsample(canonical_unitations, SizePalmerRanks);
    NullDistr(i) = ktaub([test_sample' test_sample2'],0.05,0); 
end
NullDistr = sort(NullDistr);
for i= 1:total_pseudo_replicates
    if ObsTaub < NullDistr(i)
        break
    end
end
pval1 = ((total_pseudo_replicates - i) / total_pseudo_replicates)

pval1 = pval1*total_pseudo_replicates;
plot([total_pseudo_replicates:-1:1],...
    NullDistr(:),'--')
hold on
plot([pval1,pval1],[-.2,ObsTaub],'r')
plot([0,pval1],[ObsTaub,ObsTaub],'r')
set(gca, 'FontSize',24);
%xticks([1:total_pseudo_replicates/5:total_pseudo_replicates]);
yticks([-0.2:.1:.2]);

% Create xlabel and ylabel
X = xlabel({'Sorted Test Statistic Index'});
Y = ylabel({'Sorted Test Statistic Value'});
set(X, 'FontSize', 24, 'FontName','Helvetica');
set(Y, 'FontSize', 24, 'FontName','Helvetica');
%set(gca,'xscale','log')


% gcf is some kind of magic handle to the figure; I found it in the
% saveas help page w/o much explanation.
saveas (gcf,'../PalmerFig2Fig','pdf');

hold off
