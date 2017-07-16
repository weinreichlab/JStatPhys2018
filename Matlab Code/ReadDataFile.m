% Code to read in data csv-formated files. Syntax:
%
% * zero or more comment lines, starting with #
% * one line containing the cardinalities of each locus. Code assumes
% biallelic loci, and only counts the number of cardinalities (in order to
% get the number of loci).
% * optionally a line containing the fraction of model variance due to
% experimental variance
% * one line containing phenotype names. Here parenthesis or another #
% means comment. (That allows me to include variance data if I like, but
% also use the # semantic from above.)
% * 2^loci lines of data: a series of loci 0's and 1's and then the
% phenotype(s)

function [genotypes,phenotype_count,phenotype_name,...
    phenotype_column,phenotype_experimental_error,csv] = ...
    ReadDataFile(file,directory)
    % build the filename
    filename = strcat(directory,file.name);
    % open the file and parse out the headers
    fileID = fopen(filename);
    
    % skip over the comments
    comment_lines = 0;
    while (1)
        tline = fgetl(fileID);
        if tline(1) == '#'
            comment_lines = comment_lines + 1;
        else
            break;
        end
    end
    
    % next line in file has cardinalities, which we can use to get # bits
    % (though this implementation assumes no two-digit cardinalities).
    
    % And more seriously, Walsh decomposition can't handle any cardinalites
    % larger than 2, though I don't explicitly check for that condition...
    bits = size(tline,2);
    % guaranteed to be at least one trailing comma because there has to be
    % at least one phenotype column. (To be quite precise, csv format has
    % a comma only between each column. But it makes sure that there are
    % the same number of commas on every line. In this case we know that
    % there will be at least one comma after the last cardinality value
    % because the actual data lines have one column per locus plus at least
    % one phenotype column.) We want to walk backwards from the end of the
    % line to find that first traiing comma.
    while (1)
        if tline(bits) == ','
            bits = bits - 1;
        else break
        end
    end
    % this walked us back to the last digit, which is one space too many
    bits = (bits+1)/2;
    
    % And this can be converted to the size of the dataset
    genotypes = 2^bits;
    
    % Now get the name(s) of (all) the phenotype(s).
    % This is slightly kludgy, in that phenotypes are recognized by NOT
    % having any parens in the name. In contrast, variance, SD, STDERR and
    % CI's are all required (by me!) to have parens in the column heading.
    tline = fgetl(fileID);
    
    % this little beauty loads each element of this line into a cell array
    C=textscan(tline,'%s','delimiter',',');
    % (cell array requires use of curly braces, which I don't really
    % understand. Hence the somewhat convoluted passage through z in the
    % loop...)
    C_size = size(C{1},1);
    phenotype_count = 0;
    
    
    for i = bits+1:C_size
        z = C{1}(i);
        % for some to-me stupid MatLab reason I can't just take a logical
        % or.
        if strfind(z{1},'(') 
            continue
        end
        if strfind(z{1},'#')
            continue
        end
        phenotype_count = phenotype_count + 1;
        phenotype_name(phenotype_count) = C{1}(i);
        phenotype_column(phenotype_count) = i;
        phenotype_experimental_error(phenotype_count) = 0;
    end

    fclose(fileID);
    
    % read data: skip over the comments plus cardinality and locus name lines
    csv = csvread(filename,comment_lines+2);
    
    % We optionally store experimental variance on the first line of the
    % csv
    if size(csv,1) ~= genotypes
        for i=1:phenotype_count
            phenotype_experimental_error(i) = ...
                csv(1,phenotype_column(i));
        end
        % In that case, strip off the first line
        csv = csv(2:genotypes+1,:);
        % and convert to percent
        phenotype_experimental_error = phenotype_experimental_error .* 100;
    end
end
