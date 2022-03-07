function varargout = qNorFACSimg(FACSdata, IMGdata, channelsOfInterest,varargin)

% This function is to normalize FACS data and image data using quantile
% approach (using image data as reference distribution). After normalization, 
% the distribution of relevant channel in FACSdata and IMGdata is identical. 
% Meanwhile, the rank of each entry (cell or pixel) in the relevant channel
% is maintained.
%
% Usage: [qNorFACSdata, qNorIMGdata] = qNorFACSimg(FACSdata, IMGdata,channelsOfInterest)
%        qNorFACSimg(FACSdata, IMGdata,channelsOfInterest)
%        qNorFACSimg(...,'Transfer',0 or 1)
%
% Input
% The input FACSdata should be a n*p array (n: bacterial cells; p: number
% of channels).
%
% The input IMGdata should be a m*p array (m: image pixels; p: number
% of channels).
%
% The input channelsOfInterest is a cell array, specifying p channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'.
%                        
% The optional input 'Transfer' should be either 0 
%  (when there is NO reference distribution)
%                                            or  1 (default)
% (when the IMG data serve as the reference distribution)     
%
% Output
% The optional output qNorFACSdata is a n*p array, with identical
% distribution of another optional output IMGdata (m*p) for each of the p
% channels.
%
% If the output is not assigned to any variable, the function shows the
% histogram of all channels before and after normalization (p*2 plots) the
% histogram of relevant channels after normalization (p plots); the scatter
% plot of all entries of all channels before and after normalization (p
% plots, used to check monotonicity);
%
% Dependency: quantilenormT 
% (The original form of 'quantilenormT' is Mathworks' function 
% 'quantilenorm', the modification is made for reference distribution 
% specification.)
%
% Written by Jiejue Li @ Tsinghua University
% Version 0.0.3. Created on Sep 26th, 2019. Last modified on Oct 4th, 2019.

% optional
argin = inputParser;
addParameter(argin,'Transfer',1)
parse(argin,varargin{:})
Transfer = argin.Results.Transfer;

N = size(FACSdata,1);
M = size(IMGdata,1);
NormFACS = NaN(size(FACSdata)); 
NormIMG = NaN(size(IMGdata));

for p=1:numel(channelsOfInterest)
    % split to fit the input requirement of the function 'quantilenormT'
    Splitchan = NaN(max(N,M),2);
    Splitchan(1:N,1) = FACSdata(:,p);
    Splitchan(1:M,2) = IMGdata(:,p);
    % quantilenorm
    pnormalizeddata = quantilenormT(Splitchan,'Transfer',Transfer); 
    % after normalization, output results to the new Normilized matrix 
    NormFACS(:,p) = pnormalizeddata(1:N,1);
    NormIMG(:,p) = pnormalizeddata(1:M,2);
end

% output
if nargout == 0
    close all;
    figure('Name','Histogram of all channels before and after normalization')
    for m = 1:numel(channelsOfInterest)
        subplot(2,4,m)
        h1=histogram(NormFACS(:, m), 'DisplayStyle','bar', 'FaceColor','r','FaceAlpha',0.3);
        hold on
        h2=histogram(FACSdata(:, m),'DisplayStyle','bar', 'FaceColor','g','FaceAlpha',0.3);
        h1.BinWidth = 0.01;
        h2.BinWidth = 0.01;
        legend('Normalized','Original')        
        title([channelsOfInterest{m},' FACS'])
    end
    
    for m = 1:numel(channelsOfInterest)
        subplot(2,4,m+4)
        h1=histogram(NormIMG(:, m), 'DisplayStyle','bar', 'FaceColor','r','FaceAlpha',0.3);
        hold on
        h2=histogram(IMGdata(:, m),'DisplayStyle','bar', 'FaceColor','g','FaceAlpha',0.3);
        h1.BinWidth = 0.01;
        h2.BinWidth = 0.01;
        legend('Normalized','Original')
        title([channelsOfInterest{m},' IMG'])
    end
    
    figure('Name','histogram of relevant channels AFTER normalization')
    for m = 1:numel(channelsOfInterest)
        subplot(2,2,m)
        h1=histogram(NormFACS(:, m), 'DisplayStyle','bar', 'FaceColor','b','FaceAlpha',0.3);
        hold on
        h2=histogram(NormIMG(:, m),'DisplayStyle','bar', 'FaceColor','m','FaceAlpha',0.3);
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.BinWidth = 0.01;
        h2.BinWidth = 0.01;
        title(channelsOfInterest{m})
        legend('FACS Norm','IMG Norm')        
    end
    
    figure('Name','histogram of relevant channels BEFORE normalization')
    for m = 1:numel(channelsOfInterest)
        subplot(2,2,m)
        h1=histogram(FACSdata(:, m), 'DisplayStyle','bar', 'FaceColor','b','FaceAlpha',0.3);
        hold on
        h2=histogram(IMGdata(:, m),'DisplayStyle','bar', 'FaceColor','m','FaceAlpha',0.3);
        h1.Normalization = 'probability';
        h2.Normalization = 'probability';
        h1.BinWidth = 0.01;
        h2.BinWidth = 0.01;
        legend('FACS Original','IMG Original')        
        title(channelsOfInterest{m})
    end
    
    figure('Name','Monotonicity Check')
    for m = 1:numel(channelsOfInterest)
        subplot(2,4,m)
        scatter(FACSdata(:,m),NormFACS(:,m),5,'filled')
        title([channelsOfInterest{m} ' FACS'])
        xlabel('Original values')
        ylabel('Normalized values')
    end
    
    for m = 1:numel(channelsOfInterest)
        subplot(2,4,m+4)
        scatter(IMGdata(:,m),NormIMG(:,m),5,'filled')
        title([channelsOfInterest{m} ' IMAGE'])
        xlabel('Original values')
        ylabel('Normalized values')
    end
    
else
    varargout{1} = NormFACS;
    varargout{2} = NormIMG;
end
end


