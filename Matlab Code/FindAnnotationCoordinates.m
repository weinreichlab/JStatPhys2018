% This piece of code was written to empirically find what coordinates to
% use on my annotations of the Experimental Error line in these figures.
%
% DMW June 15, 2017
%
% It seems that 10% of the figure is devoted to the margins, so to place
% the label at position x means specifying it be at 0.8x+0.1.
%
% DMW June 23, 2017
% I want the aspect ratio to be the golden mean (thanks Stephen Wolfram!)
% which requires recalculating this. I didn't bother to fully suss it out
% this time. Instead, the first block has a function that was tweaked
% empirically to work.
clear all
close all

length = .75;
height = 0.0651;
max = 10;
for i=0:max
    plot([0,1],[i/max,i/max],'--');
    pbaspect([1.618 1 1]);
    hold on
    % And label it. 
    t=annotation('textbox','string',sprintf('Experimental variance %f',i/max));
    t.FontSize =12; 
    t.LineStyle = 'none';
    t.Position = [0.15, .18+i*1.025/(1.618*max),length,height];
end
hold off

% Test those results but then needed to tweak to 0.09+0.815*x
% figure (2)
% for i=0:max
%     plot([0,1],[i/max,i/max],'--');
%     pbaspect([1.618 1 1]);
%     hold on
%     % And label it. 
%     t=annotation('textbox','string',sprintf('Experimental variance %f',i/max));
%     t.FontSize =12; 
%     t.LineStyle = 'none';
%     t.Position = [0.15, 0.09+0.815*i/max,length,height];
% end
% hold off

