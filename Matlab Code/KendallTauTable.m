% Code to extract empirical magnitude-sorted Walsh coefficient ranks,
% compute Kendall's tau_b, and then compare to a null-hypothesis
% distribution of values generated for random Walsh coefficient ranks.
%
% Derived from Jacob Jaffe's code, downloaded from GitHub June 7, 2017.
%
% DMW June 16, 2017

clear all;

files = dir('../Datasets/*.csv');
directory = '../Datasets/';

filenames = [];
phenotype_names = [];
number_of_loci = [];
number_of_peaks = [];
number_of_Walsh_coefficients = [];
tau_b = [];
pvalue_1 = [];

% Building these null distributions is fairly slow so I wrote a piece of
% code (BuildNullTau_bDistr.m) to do it once and save the results in a 
% file.
MAX_DATASET_SIZE = 6;
load('./NullTau_bDistr.mat');
total_pseudo_replicates = size(taub_null_distr,1);

% loop through all the data files
for file = files'
    % This abstracts out reading the datafiles for use here and in
    % the code that draws the Residual Variance figures
    [genotypes,phenotype_count,phenotype_name,...
        phenotype_column,phenotype_experimental_error,data] = ...
        ReadDataFile(file,directory);
    
    % Capture this for the output file
    filename = string(file.name);
    
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
%   does.
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
        
        % We previously will have built the null distribution for datasets
        % no larger than MAX_DATASET_SIZE. So we need to make sure we don't
        % blow that limit.
        if number_of_loci(end) > MAX_DATASET_SIZE
            disp('ERROR: Dataset too large');
            break;
        end
        
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
        for i=1:genotypes
            if Residual_Variance(i) < ...
                phenotype_experimental_error(phenotype_index)
                i = i - 1;
                break;
            end
        end
        truncated_sorted_empirical_unitations = ...
            sorted_empirical_unitations(1:i);
        truncated_canonical_unitations = ...
            canonical_unitations(1:i);
        
        % Compute tau_b for this phenotype compared to canonical (naive)
        % unitation order. Record it and the number of coefficients used in
        % the computation. 

        tau_b(end+1) = ...
            ktaub([truncated_sorted_empirical_unitations ...
            truncated_canonical_unitations],0.05,0);
        number_of_Walsh_coefficients(end+1) = i;

        % Now find the p-value for this dataset.
        if number_of_Walsh_coefficients(end) > 1
            for i= 1:total_pseudo_replicates
                if tau_b(end) < taub_null_distr(i,number_of_Walsh_coefficients(end))
                    break
                end
            end
        else
            % If there aren't at least 2 significantly non-zero Walsh
            % coefficients then the p-value is 1.
            i = 0;
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
writetable(Table,'../TauTable.csv') ;