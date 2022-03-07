function normalizedVals = quantilenormT(values,varargin)
% QUANTILENORM performs quantile normalization over multiple arrays
%
%   NORMDATA = QUANTILENORM(DATA), where the columns of DATA correspond to
%   separate chips, normalizes the distributions of the values in each
%   column. Note that if DATA contains NaN values, then NORMDATA will also
%   contain NaNs at the corresponding positions.
%    
%   optional input:
%   quantilenormT(...,'Transfer',0 or 1)
%    The optional input 'Transfer' should be either 0 
%   (when there is NO reference distribution)
%                                            or     1 (default)
%   (when the SECOND column in data serve as the reference distribution)       
%
%   Examples:
%       load yeastdata
%       normYeastValues = quantilenorm(yeastvalues,'display',1);
%
%   See also AFFYGCRMA, AFFYINVARSETNORM, AFFYPREPROCESSDEMO, AFFYRMA,
%   GCRMA, GCRMABACKADJ, MAINVARSETNORM, MALOWESS, MANORM, RMABACKADJ,
%   RMASUMMARY.

% Reference:
% Probe Level Quantile Normalization of High Density Oligonucleotide Array
% Data. Bolstad, B. http://stat-www.berkeley.edu/~bolstad/stuff/qnorm.pdf
% Copyright 2004-2008 The MathWorks, Inc.
% Last modified on Sep 26th, 2019. by Jiejue Li @ Tsinghua University

argin = inputParser;
addParameter(argin,'Transfer',1)
parse(argin,varargin{:})
Transfer = argin.Results.Transfer;


% bioinfochecknargin(nargin,1,mfilename);
% % set defaults
tiedFlag = true;
dispFlag = false;
avgFcn = @nanmean;
% % deal with the various inputs
% if nargin > 1
%     if rem(nargin,2) == 0
%         error(message('bioinfo:quantilenorm:IncorrectNumberOfArguments', mfilename));
%     end
%     okargs = {'ignoreties','display','median'};
%     for j=1:2:nargin-2
%         pname = varargin{j};
%         pval = varargin{j+1};
%         k = find(strncmpi(pname,okargs,numel(pname)));
%         if isempty(k)
%             error(message('bioinfo:quantilenorm:UnknownParameterName', pname));
%         elseif length(k)>1
%             error(message('bioinfo:quantilenorm:AmbiguousParameterName', pname));
%         else
%             switch(k)
%                 case 1  % tied flag
%                     % quick and dirty method when there are not many tied
%                     % values. This may be redundant as tiedrank is pretty
%                     % efficient.
%                     tiedFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
%                 case 2 % display flag
%                     dispFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
%                 case 3 % median flag
%                     medianFlag = bioinfoprivate.opttf(pval,okargs{k},mfilename);
%                     if medianFlag
%                         avgFcn = @nanmedian;
%                     end
%             end
%         end
%     end
% end

% allocate some space for the normalized values
normalizedVals = values;
valSize = size(values);
rankedVals = NaN(valSize);
% find nans
nanvals = isnan(values);
numNans = sum(nanvals);
ndx = ones(valSize);
N = valSize(1);
% create space for output
if tiedFlag
    rr = cell(size(values,2),1);
end
% for each column we want to ordered values and the ranks with ties
for col = 1:valSize(2)
    [sortedVals,ndx(:,col)] = sort(values(:,col));
    if(tiedFlag)
        rr{col} = sort(tiedrank(values(~nanvals(:,col),col)));
    end
    M = N-numNans(col);
    % interpolate over the non-NaN data to get ranked data for all points
    rankedVals(:,col) = interp1(1:(N-1)/(M-1):N,sortedVals(1:M),1:N,'spline');
end
% take the mean of the ranked values
if Transfer==0
mean_vals = feval(avgFcn,rankedVals,2);
end
if Transfer==1
mean_vals = rankedVals(:,2);    %hold the second row data
end
% display the estimated distributions of the columns
if dispFlag
    numBuckets = min(valSize(1)/10,30);
    x = linspace(min(mean_vals), max(mean_vals),numBuckets);
    n = histc(mean_vals,x);
    nAll = repmat(n,1,valSize(2));
    for count = 1:valSize(2)
        x = linspace(min(mean_vals), max(mean_vals),numBuckets);
        nAll(:,count) = histc(values(:,count),x);
    end
    % end point from histc will always be 1 for n and probably 0 for nAll
    % so don't show this in case it makes things look ugly.
    plot(x(1:end-1),nAll(1:end-1,:)); hold on;
    plot(x(1:end-1),n(1:end-1,:),'k','lineWidth',3);
    legendString = strread(sprintf('Distribution %d\n',(1:valSize(2))'),...
        '%s','delimiter','\n');
    legendString{end+1} = 'Normalized Distribution';
    legend(legendString);
    hold off;
end

% Extract the values from the normalized distribution
for col = 1:size(values,2)
    M = N-numNans(col);
    if tiedFlag
        normalizedVals(ndx(1:M,col),col) = interp1(1:N,mean_vals,1+((N-1)*(rr{col}-1)/(M-1)));
    else
        normalizedVals(ndx(1:M,col),col) = interp1(1:N,mean_vals,1:(N-1)/(M-1):N,'spline');
        
    end
end

