% Code to extract empirical magnitude-sorted Walsh coefficient ranks,
% compute Kendall's tau_b, and then compare to a null-hypothesis
% distribution of values generated for random Walsh coefficient ranks.
%
% Derived from Jacob Jaffe's code, downloaded from GitHub June 7, 2017.
%
% DMW June 16, 2017
%
% DMW July 14, 2017: just realized that at present the null distribution 
% for x significantly-residual-variance-reducing Walsh coefficients assumes
% that those epistatic ranks are sampled from those seen for exactly
% ceil(log2(x)) loci. For example, if I have 6 rank orders to compare to
% the canonical expectation, I assume that the canonicals are the first six
% of those defined by ceil(log2(6)) = 3 bits, i.e., 0 1 1 1 2 2 (2 3),
% where the (2 3) are those lopped off.
%
% But if the 6 significantly-residual-variance reducing Walsh coefficients
% are sampled from those seen for, say 5 loci, then the first six
% canonicals would be 0 1 1 1 1 1 1, and permuting those would give quite a
% different null distribution of tau_b values against which to assess
% significance.
%
% So I /could/ pre-fabricate a big table of tau_b for permuted values of x
% significantly-residual-variance-reducing coefficients for all possible
% values of y loci, but that's crazy! Instead, I'll build the distr's for
% just those values of x and y that turn up in the data right here. But
% each that I encounter I also save in a file for future reference.

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
        
        % get the phenotype vector
        phenotype = data(:,phenotype_column(phenotype_index));
        % Get the number of peaks on this landscape
        number_of_peaks(end+1) = CountPeaks(genotypes,phenotype);
        % Get the Walsh coefficients for this phenotype
        Walsh = hadamard(genotypes) * phenotype/genotypes;
    
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
        
        % Compute tau_b for this phenotype compared to canonical (naive)
        % unitation order. Record it and the number of coefficients used in
        % the computation. 

        tau_b(end+1) = ...
            ktaub([truncated_sorted_empirical_unitations ...
            truncated_canonical_unitations],0.05,0);
        number_of_Walsh_coefficients(end+1) = sig_nonzero_Walsh;

        % Now find the p-value for this dataset. First build the filename
        NullTaubFilename = ...
            sprintf('./NullTau_bDistrs/%dLoci%dSigWalsh.mat',loci,...
            sig_nonzero_Walsh);
        % Then see if it already exists
        if exist(NullTaubFilename)
            % if so, use it
            load(NullTaubFilename)
            total_pseudo_replicates = size(taub_null_distr,1);
        else
            % if not, build it: permute the truncated canonical unitations.
            total_pseudo_replicates = default_pseudo_replicates;
            taub_null_distr = zeros( total_pseudo_replicates, 1);
            for i = 1:total_pseudo_replicates
                test_sample = randsample(truncated_canonical_unitations, ...
                    sig_nonzero_Walsh);
                test_sample2 = randsample(truncated_canonical_unitations, ...
                    sig_nonzero_Walsh);
                taub_null_distr(i) = ktaub([test_sample test_sample2],...
                    0.05,0);
            end
            taub_null_distr = sort(taub_null_distr);
            % And save it for next time
            save(NullTaubFilename,'taub_null_distr');
        end
        
        % Finally, find where our observation fits among the permutations
        for i= 1:total_pseudo_replicates
            if tau_b(end) < taub_null_distr(i)
                break
            end
        end
        
        % Since our naive hypothesis is that the sorted_empirical_untation
        % should be perfectly correlated with canonical_unitations, we want
        % to know what fraction of the pseudo-random replicates have a
        % value larger than the empirical one. 
        pval1 = ((total_pseudo_replicates - i) / total_pseudo_replicates);
        pvalue_1(end+1) = pval1;
    end
end

% Now build a table containing these results and save as a .csv file.
% First transpose all the row vectors into columns.
filenames = filenames';
phenotype_names = phenotype_names';
number_of_loci = number_of_loci';
number_of_peaks = number_of_peaks';
number_of_Walsh_coefficients = number_of_Walsh_coefficients';
tau_b = tau_b';
pvalue_1 = pvalue_1';

% Build the Bonferroni threshold p-values, which depends on the total
% number of tests. 
table_size = size(tau_b,1);
Bonferroni_factor = zeros(table_size,1);
for i = 1:table_size
    Bonferroni_factor(i) = 1/(table_size-i+1);
end

Bonferroni_05 = 0.05*Bonferroni_factor;
Bonferroni_01 = 0.01*Bonferroni_factor;
Bonferroni_001 = 0.001*Bonferroni_factor;

% Assign corrected p-values
Corrected_p = strings(table_size,1);
sort_pvalue_1 = sort(pvalue_1);
for i = 1:table_size
    if sort_pvalue_1(i) < Bonferroni_001
        Corrected_p(i) = '***';
    elseif sort_pvalue_1(i) < Bonferroni_01
        Corrected_p(i) = '**';
    elseif sort_pvalue_1(i) < Bonferroni_05
        Corrected_p(i) = '*';
    end
end

% makes and sorts table by pvalues in ascending order
Table = table(filenames, phenotype_names, number_of_loci, ...
    number_of_peaks,number_of_Walsh_coefficients, tau_b,  pvalue_1);
Table = sortrows(Table, 'pvalue_1');

% paste in the Bonferroni critical values and the observed significance
% asterisks. MODIFIED TO ONLY SAVE CORRECTED P-VALUES.
% Table = [Table, table(Bonferroni_05),...
%     table(Bonferroni_01), table(Corrected_p)];
Table = [Table, table(Corrected_p)];

% will write Table of statistics to a csv file
if NK
    writetable(Table,'../NKTauTable.csv') ;
else
    writetable(Table,'../TauTable.csv');
end