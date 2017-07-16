% Code to draw data figs for JStatPhys paper. Derived from Jacob Jaffe's 
% implementation (Simplified_Regression_Main_with_plots.m), downloaded 
% from github June 7, 2017. 
%
% Change aspect ratio to golden mean (which also required re-jiggering the
% math for placing legend and text on Experimental Variance line) and
% increase font size to allow for shrinking. I reduced Palmer2015 by 50% in
% Illustrator when I made Figure 1.
%
% DMW June 23, 2017

clear all
close all

% This flag determines whether the data are log-transformed before analysis
% In the CODEG paper we log-transformed drug resistance and growth rate
% datasets, but here I've put the raw data in the dataset files and will
% run this analysis twice.
%
% One slightly kludgy consequence of this implementation is that if 
% all the datasets in a data file have at least one zero then the code 
% produces a file but there's nothing plotted. At present, that's Brigdham
% and Palmer.
%
% OBSOLETED JUNE 13, 2017: INDIVIDUAL DATA FILES ARE NOW LOG-TRANSFORMED OR
% NOT AS PER CODEG 2013 CONVENTION.

LOG_TRANSFORM_DATA = 0;

% dynamic range of y-axis.
y_min = 10^-4;
y_max = 10^2;
% load the color codes for plotting
colors = {'k','b','m','r','g','y','c'};
num_colors = size(colors,2);

% Do we want data from the simulated control or the biology?
NK = 0;

if NK
    files = dir('../NK landscapes/*.csv');
    directory = '../NK landscapes/';
else
    files = dir('../Datasets/*.csv');
    directory = '../Datasets/';
end

% I changed something in Jacob's code, after which it started overwriting
% each file's results onto the same panel. So I now use an index and call
% figure(figure_number) before drawing each figure.
figure_number = 1;

%loop through all the *.csv datafiles in the directory
for file = files'
    % This abstracts out reading the datafiles for use here and in
    % the code that does the analysis of Kendall's tau_b.
    [genotypes,phenotype_count,phenotype_name,...
        phenotype_column,phenotype_experimental_error,csv] = ...
        ReadDataFile(file,directory);
       
    % open the window
    figure(figure_number)
    % Create xlabel and ylabel
    X = xlabel({'Number of Parameters'});
    Y = ylabel({'Residual Variance (%)'});
    set(X, 'FontSize', 24, 'FontName','Helvetica');
    set(Y, 'FontSize', 24, 'FontName','Helvetica');
    
    % We'll use this to set the lower bound on the y axis.
    % Though in fact we've standardized y-axis dynamic range.
%     ymax = 100;
%     ymin = ymax;
    
    for phenotype_index = 1:phenotype_count
        phenotype = csv(:,phenotype_column(phenotype_index));
        
        % This block of code log-transforms the data before going to work
        % First thing is to protect from any phenotype values == 0, which 
        % would cause the log-transformation to blow up. In that case, skip 
        % this phenotype.
        if LOG_TRANSFORM_DATA == 1
%            if prod(csv(:,phenotype_column(phenotype_index))) == 0
            if all(phenotype > 0) == 0
                t = annotation('textbox');
                t.String = ...
                    'At least one zero value so unable to log-transform';
                continue
            end
            phenotype = log(phenotype);
        end
    
        % function returns values listed below from simplified regression
        [Rsquared, unitation] = ...
            GetResidualVariance(phenotype);
      
        % For graphing residual variance. Second line a bit of a kludge
        % but for some reason, for some datasets fitlm() return a 1 for 
        % the first R^2, instead of 0.
        Residual_Variance = 1-Rsquared;
        Residual_Variance(1) = 1;
        % Convert to percent
        Residual_Variance = Residual_Variance * 100;
    
        % We don't want to try to plot any values below the y-axis...
        for i=1:genotypes
            if Residual_Variance(i) < y_min
                break
            end
        end
        % i points to the first too-small value, so back it up by one.
        i = i-1;
        % This beauty does the plotting, labeling with each point's unitation
        text([1:i]',Residual_Variance(1:i),num2str(unitation(1:i)),...
            'color',colors{mod(phenotype_index-1,num_colors)+1},...
            'FontSize',12);
        hold on
        pbaspect([1.618 1 1]);
        set(gca, 'FontSize',18);
        
        % This horror buids the legend
        xbegin = 0.17;
        ybegin = 0.22;
        % these two are fairy unclear to me but seem by trial to error to
        % work
        length = .9;
        height = 0.0651;
        t = annotation('textbox');
        t.String = phenotype_name{phenotype_index};
        t.Color = colors{mod(phenotype_index-1,num_colors)+1};
        t.LineStyle='none';
        t.FontSize=18;
        t.Position=[xbegin,ybegin + ...
            (phenotype_count - phenotype_index)*0.05,length,height];
        
        % Set y axis bounds. 
        % First, find the smallest non-zero value
%         for i=genotypes:-1:1
%             if Residual_Variance(i) > 0
%                 break
%             end
%         end
        % and if this phenotype's minimum is the smallest seen so far, use
        % it
%         if Residual_Variance(i) < ymin
%             % Then round down to the nearest power of 10
%             y_lowerbound = 10^floor(log10(Residual_Variance(i)));
%         end
        
        % THOUGH ACTUALLY I'VE DECIDED TO STANDARDIZE TO y_min = 10^-4 and 
        % y_max = 10^2. See above for variable declarations and
        % initializations
       
        % Now plot the proportion of the model variance that's represented
        % by the experimental variance if known
        if phenotype_experimental_error(phenotype_index) > 0
            plot([0,genotypes],[phenotype_experimental_error(phenotype_index),...
                phenotype_experimental_error(phenotype_index)],'--','color',...
                colors{mod(phenotype_index-1,num_colors)+1})
            % And label it. 
            t=annotation('textbox','string','Experimental variance');
            t.Color = colors{mod(phenotype_index-1,num_colors)+1};
            t.LineStyle = 'none';
            t.FontSize =18; 
            % This business is a bit ugly! the coordinates passed represent
            % the proportions across the figure: [0,] = lower left corner
            % and [1,1] represents upper right. To convert to coordinates
            % within the actual field of the plot I wrote
            % FindAnnotationCoordinates.m, which shows that the margins of
            % the plots account for the first and last 10% of the plot.
            % That would suggest a conversion formula of 0.1+0.8*y to get
            % the y I want. But add to that some unknown padding around the
            % text and I conclude that the correct formula is 0.09+0.815*y.
            %
            % Things are then made slightly more complicated by the fact
            % that we log-transform the y-axis. So finally we need
            y = 0.09 + 0.815* ...
                log10(phenotype_experimental_error(phenotype_index)/y_min)/...
                log10(y_max/y_min);
            % All the tidy foregoing logic went out the window when I
            % changed to the Mathematican golden aspect ratio... What
            % follows is purely empirical.
            y = 0.19 + 1.025/1.618* ...
                log10(phenotype_experimental_error(phenotype_index)/y_min)/...
                log10(y_max/y_min);
            t.Position = [0.17,y,length,height];
        end
    end
    
    % do scaling
    xlim( [0,genotypes]);
    xticks([0:genotypes/8:genotypes]);
    % I'd like to standardize this for ease of comparision across panels,
    % so y_lowerbound actually isn't used for scaling the y axis.
%     if y_lowerbound < y_min
%         y_lowerbound = y_min;
%     end
    
    ylim( [y_min,y_max]);
    yticks(10.^[log10(y_min):log10(y_max)]);
    set(gca,'yscale','log')
    New_YTickLabel = get(gca,'ytick');
    set(gca,'YTickLabel',New_YTickLabel);

    % This is probably pretty bad style: overwriting a system-loaded
    % variable like file.name. But I'm lazy!
    filename_length = size(file.name,2);
    file.name(filename_length-2:filename_length) = 'pdf';
    if LOG_TRANSFORM_DATA == 0
        output_filename = sprintf( '../Residual Variance Figs/%s',file.name);
    else
        output_filename = sprintf( '../Figure 1 ln(data)/%s',file.name);
    end
    
    if ~exist('../Residual Variance Figs')
        mkdir('../Residual Variance Figs');
    end
    
    % gcf is some kind of magic handle to the figure; I found it in the
    % saveas help page w/o much explanation.
    saveas (gcf,output_filename,'pdf');

    hold off
    figure_number = figure_number+1;
end
