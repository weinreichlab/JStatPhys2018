% this function is a simplified version of a stepwise linear regression. It
% takes advantage of the fact that Walsh coefficients with greater absolute
% magnitued will have more explanatory power in a regression model.
%
% Code almost verbatim from Jacob Jaffe's implementation, downloaded from
% github June 7, 2017. Two difference:now there's only one input argument 
% (instead I compute Walsh coefficients here) and I dropped his aicc code 
% (which wasn't applicable to these data for a reason I don't quite recall).
%
% DMW June 13, 2017

function [r_squared, sorted_unitations] = ...
    GetResidualVariance(phenotype)
    l = length(phenotype);
    H = hadamard(l);
    Walsh = H * phenotype/l;
    %
    % returns the rank of each walsh coefficient
    %
    % unique() sounds like a bad choice: what if two (or more) Walsh
    % coefficients are exactly equal. But in fact the third vector returned
    % holds the indices into the argument vector, sorted from smallest to
    % largest BEFORE DUPS ARE ELIMINATED.
    [~, ~, Walsh_ranking] = unique(-abs(Walsh));

    unitations = zeros(l, 1); % records epistatic order 
    Walsh_number = zeros(l, 1); % records indices of walsh coefficients, used in sorting
    
    
    % simple loop to get epistatic orders of each walsh coefficient
    for i= 0:l-1
        bin = dec2bin(i); %converts index to binary
        b = str2num(bin(:)); %converts binary strings to separate 0's and 1's
        unitations(i+1) = sum(b); %sums the ones to get Walsh coeffient's epistatic order
        Walsh_number(i+1) = i+1; %records the corresponing walsh number (not in bitstrings)
    end
    
    % matrix used to sort Walsh coefficients, Walsh numbers and unitations
    % based on ranking by absolute maginuted
    sort_matrix = [Walsh_ranking, Walsh, Walsh_number, unitations];
    sort_matrix = sortrows(sort_matrix, 1);

    % is subset of all walsh coefficients, iteration of regression the next 
    %walsh coefficient of greatest absolute magnitude is added to this subset
    test_walsh_coefficients = zeros(l,1); 
    
    r_squared = zeros(l,1); % records r_squared of each model of increasing complexity
    aicc = zeros(l,1); % records aic of each model of increasing complexity

    % Iterative regression model
    % Will add the next highest Walsh coefficient by magnitude to the
    % regression.
    for j = 1:l
        
        Walsh_coefficient_index =  sort_matrix(j, 3); % gets index
        Walsh_coefficient = sort_matrix(j, 2); % gets Walsh coefficient
        
        % places walsh coefficient into the right index of the test Walsh
        % coefficient array used in the regression
        test_walsh_coefficients(Walsh_coefficient_index) = Walsh_coefficient;
        
        backtransformed_fitness = H * test_walsh_coefficients;
        
        % performs linear regression
        model = fitlm(phenotype,backtransformed_fitness); 
        
        % records r squared term
        r_squared(j) = model.Rsquared.Ordinary;
        
    end
    
    % This is the other key thing we return.
    sorted_unitations = sort_matrix(:,4);
    
    % other thing that may be useful. Not currently returned to caller...
    sorted_Walsh_number = sort_matrix(:,3);
end

