% Code to produce Table 3 and its ilk: the average per-Walsh term
% explanatory power as a function of order. It's the same loop as
% BuildKendallTauTable.m but instead of computing tau_b and assessing
% significance, it does this other thing.
%
% DMW Jan 10, 2018

clear all;

% Do we want data from the simulated control or the biology?
NK = 0;

if NK
    files = dir('../NK landscapes/*.csv');
    directory = '../NK landscapes/';
else
    files = dir('../Datasets/*.csv');
    directory = '../Datasets/';
end

filenames = [];
phenotype_names = [];
number_of_loci = [];
number_of_peaks = [];
number_of_Walsh_coefficients = [];
tau_b = [];
pvalue_1 = [];

default_pseudo_replicates = 1e5;

% loop through all the data files
for file = files'
    % This abstracts out reading the datafiles for use here and in
    % the code that draws the Residual Variance figures
    [genotypes,phenotype_count,phenotype_name,...
        phenotype_column,phenotype_experimental_error,data] = ...
        ReadDataFile(file,directory);
    
    % Capture this for the output file
    filename = string(file.name);
    disp(sprintf('%s',filename));
    
    % Creates unitation list in canonical Walsh coefficient ordering for
    % this sized dataset.
    
%   This was Jacob's code, which gives the unitations sorted by canonical
%   Walsh indices.
%    canonical_unitations = zeros(genotypes, 1);
%     for i= 0:genotypes-1
%         % Converts the decimal index to a string of binary digits
%         bin = dec2bin(i); 
%         % The (:) operator turns a string into an array of individual
%         % characters. Each element in that array is then converted to an
%         % array of (decimal) digits...
%         b = str2num(bin(:)); 
%         % ...which can be summed to get the coefficient's rank.
%         canonical_unitations(i+1) = sum(b); 
%     end
%
%   But what we want is these same unitations, sorted by expected
%   magnitudes of associated Walsh coefficients. That's what this code
%   does instead.
    canonical_unitations = [];
    loci = log2(genotypes);
    for i=0:loci
        for j=1:nchoosek(loci,i)
            canonical_unitations(end+1) = i;
        end
    end
    canonical_unitations = canonical_unitations';

    % Cycle throught the >=1 phenotypes in this datafile
    for phenotype_index = 1:phenotype_count
        % Store this for the output file. I think this use of cat() is
        % synonymous with the use of <array>(end) used elsewhere in this
        % code.
        filenames = cat(2,filenames, filename);
        phenotype_names = ...
            cat(2,phenotype_names,phenotype_name(phenotype_index));
        fprintf('%s\n',phenotype_name{phenotype_index});
        % get the phenotype vector
        phenotype = data(:,phenotype_column(phenotype_index));
        
        % Get the number of peaks on this landscape
        number_of_peaks(end+1) = CountPeaks(genotypes,phenotype);
    
        % At the end, number_of_loci() will be a vector carrying the number 
        % of  loci in each dataset. <array>(end) is the last entry in the 
        % array, so <array>(end+1) adds one more element to the end.
        number_of_loci(end+1) = loci;
           
        % This subroutine gets the sorted unitations and Rsquared's in the
        % data. We need the latter in those datasets for which we have 
        % experimental variance information so we know where to chop 
        % off the two column vectors being compared using tau_b.
        [Rsquared, sorted_empirical_unitations] = ...
            GetResidualVariance(phenotype);
        Residual_Variance = 1-Rsquared;
        Residual_Variance(1) = 1;
        % Convert to percent because ReadDataFile() converts experimental
        % error to percent.
        Residual_Variance = Residual_Variance * 100;
            
        % If 1 - RSquared's associated with any Walsh coefficients are
        % smaller than the experimental variance, don't include them in
        % computing Kendall's tau_b. (Note that ReadDataFile() loads the
        % experimental error term with zero if the datafile doesn't have a
        % value, so this code works equally for datasets w/ and w/o error
        % info: in the latter case i = genotypes at the end of the loop.)
        for sig_nonzero_Walsh=1:genotypes
            if Residual_Variance(sig_nonzero_Walsh) <= ...
                phenotype_experimental_error(phenotype_index)
                sig_nonzero_Walsh = sig_nonzero_Walsh - 1;
                break;
            end
        end
        % We want to know how many coefficients does it take to push the
        % residual variance BELOW the threshold.
        if sig_nonzero_Walsh < genotypes
            sig_nonzero_Walsh = sig_nonzero_Walsh + 1;
        end
        
        truncated_sorted_empirical_unitations = ...
            sorted_empirical_unitations(1:sig_nonzero_Walsh);
        truncated_canonical_unitations = ...
            canonical_unitations(1:sig_nonzero_Walsh);
        
        % Following suggestion of Roy Kishony and Adam Palmer
        % I compute the average reduction in variance explained per
        % epistatic term, as a function of order.
        
        counts = zeros(loci+1,1);
        VarianceReduction = zeros(loci+1,1);
        
        % Skip the zeroth-order term, which has no explanatory power.
        % And process the sorted_empirical_unitations from 2 to genotype
        % in order that the total explanatory power sum to 100%.
        %
        % In the end I ran this loop only to sig_nonzero_Walsh, meaning
        % cumulative variance explained across terms may not sum to 1. In
        % practice, of the datasets included in Table 3, only Palmer et al
        % *has* experimental variance, and their noise is so low that the
        % sum of explained variance is 99.687141 ≈ 100!.
        for i=2:sig_nonzero_Walsh
            counts(sorted_empirical_unitations(i)+1) = ...
                counts(sorted_empirical_unitations(i)+1) + 1;
            Delta = Residual_Variance(i-1) - Residual_Variance(i);
            VarianceReduction(sorted_empirical_unitations(i)+1) = ...
                VarianceReduction(sorted_empirical_unitations(i)+1) + ...
                Delta;
        end
        
        % Build table containing these results and save as .csv file.
        % Too much work to define a filename; let's just watch the data go
        % by.
        PerOrder = VarianceReduction./counts;
        Orders = (0:loci);
        [~, ~, Rank] = unique(-PerOrder);
        table(Orders',counts,VarianceReduction,PerOrder,Rank)
        fprintf('Total explanation %f\n',sum(VarianceReduction));
        fprintf('Largest ÷ 2nd-largest %f\n\n', ...
            max(PerOrder(2:end))/max(PerOrder(PerOrder<max(PerOrder(2:end)))));
    end
end
