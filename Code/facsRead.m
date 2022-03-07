function gated_cell = facsRead(FACSdata, varargin)

% This function is to read in a .fcs FACS raw data as a table
%
% Usage: gated_cell = facsRead(FACSdata)
%        gated_cell = facsRead(FACSdata, 'FSCthr',100, 'SSCthr',200)
%
% The input FACSdata should be a standard .fcs file.
%
% The optional input FSCthr and SSCthr determine the lower limit in FACS
% data to identify an event as a bacterial cell
% The default value is 100 and 200 for FSCthr and SSCthr, respectively.
%
% Written by Tianmin Wang @ Tsinghua University
% Version 0.0.1. Created on Jul 27, 2019. Last modified on Sep 26, 2019.

argin = inputParser;
addParameter(argin,'FSCthr',100)
addParameter(argin,'SSCthr',200)
parse(argin,varargin{:})
FSCthr = argin.Results.FSCthr;
SSCthr = argin.Results.SSCthr;

% struct storing common facs names
facs_mapping = struct;
facs_mapping.FSC_A = 'FSC';
facs_mapping.FSC_H = 'FSCH';
facs_mapping.FSC_W = 'FSCW';
facs_mapping.SSC_A = 'SSC';
facs_mapping.SSC_H = 'SSCH';
facs_mapping.SSC_W = 'SSCW';
facs_mapping.FITC_A = 'GFP';
facs_mapping.PE_Texas_Red_A = 'mCherry';
facs_mapping.Pacific_Blue_A = 'CFP';
facs_mapping.PE_Cy5_A = 'Cy5';
facs_mapping.APC_A = 'APC';
facs_mapping.Time = 'Time';

% read in biofilm facs data
[fcsdat, fcshdr] = fca_readfcs(FACSdata);
channels = {fcshdr.par.name};
for c = 1:length(channels) % make the name legal as column names in table calss
    channel = regexprep(channels{c},{' ','-'},{'_','_'});
    channels{c} = facs_mapping.(channel);
end
fcs_table = array2table(fcsdat, 'VariableNames',channels); % convert to table
fsc_col = table2array(fcs_table(:, 'FSC'));
ssc_col = table2array(fcs_table(:, 'SSC'));
gated_cell = fcs_table(  (fsc_col > FSCthr) & (ssc_col > SSCthr), : );