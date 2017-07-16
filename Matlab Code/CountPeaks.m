% Code to count the number of maxima on a given fitness landscape.
%
% DMW June 27, 2017

function [peaks] = CountPeaks(genotypes,phenotype)

loci = log2(genotypes);
peaks = 0;
for i=0:genotypes-1
    % this is the flag
    peak = 1;
    % Convert focal genotype to binary
    bin = dec2bin(i,loci);
    % Build every single-mutant neighbor
    for j=1:loci
        tmpbin = bin;
        if bin(j) == '1'
            tmpbin(j) = '0';
        else
            tmpbin(j) = '1';
        end
        % If focal genotype smaller than this neighbor then it's
        % not a peak. (Remember that arrays are 1-based in Matlab)
%                 display(sprintf('w(g(%d))=%f has neighbor w(g(%d))=%f',...
%                     i,phenotype(i+1),bin2dec(tmpbin),...
%                     phenotype(bin2dec(tmpbin)+1)));
        if phenotype(i+1) < phenotype(bin2dec(tmpbin)+1)
            peak = 0;
        end
    end
    % If focal was smaller than any of its neighbors then this flag
    % will have been shut off
    if peak
        peaks = peaks + 1;
%                 display(sprintf('%d is a peak',i));
    else
%                display(sprintf('%d is not a peak',i));
    end
end