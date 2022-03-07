function xmlGates = xmlGateExtract(xmlStruct, required_FACS_channels, clusterNum)

% This function is used to extract all the gate information from a xml
% struct, thus bypassing overwhelming redundant information in xml struct>
% Meanwhile, this script is also used to check whether xmlStruct is
% consistent with the required_FACS_channels and clusterNum.
% 
% Usage: xmlGates = xmlGateExtract(xmlStruct, required_FACS_channels, clusterNum)
%
% Input
% The input xmlStruct is a struct converted from a .xml file from FACS
% machine, executed by xml2struct.
%
% The input required_FACS_channels is a cell array, specifying the channels of
% interest. Accepted channels include 'Pacific_Blue_A', 'FITC_A', 
% 'PE_Texas_Red_A' and 'APC_A'.
%
% The input clusterNum is an integer, specifying the total number of cell
% subpopulations to be sorted.
%
% Output
% The output xmlGates is a struct storing the essential information about
% each gate in xmlStruct, except for 'All Events' gate.
% xmlGates.name -> xparam(str), yparam(str), xlog(logical), ylog(logical), used.
%
% Written by Ping Shen and Tianmin Wang @ Tsinghua University
% Version 0.0.1. Created on Oct 19, 2019. Last modified on Oct 19, 2019.

% extract the essential gate information from .xml file
xmlGates = struct;
gateNum = numel(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate);
FSC_SSC_gate = 'P1';
for m = 1:gateNum
    gateName = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.name.Text;
    if ~strcmp(gateName, 'All Events') % normal polygon gates
        gateType = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.region.Attributes.type;
        
        if ~strcmp('POLYGON_REGION', gateType) % gate is not polygon
            error(['Gate ', gateName, ' is not polygon!'])
        end
        
        polygonGate_points = numel(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.region.points.point);
        x_param = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.region.Attributes.xparm;
        y_param = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.region.Attributes.yparm;
        
        % check hierarchy
        if strcmp(x_param, 'FSC-A') && strcmp(y_param, 'SSC-A') % FSC-SSC gate should have 'All Events' gate as its parent
            parent_gate = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.parent.Text;
            if ~strcmp(parent_gate, 'All Events')
                
                error('FSC-SSC gate should use all events as its parent!')
            end
            FSC_SSC_gate = [parent_gate, '\', gateName];
        else % Fluorescent gate should have FSC-SSC gate as its parent and is 4-points polygon
            if polygonGate_points ~= 4
                error(['Polygon gate ', gateName, ' has more than or less than 4 points!'])
            end
            parent_gate = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.parent.Text;
            if ~strcmp(parent_gate, FSC_SSC_gate)
                error(['Fluorescent gate ', gateName, ' should use FSC-SSC gate as its parent!'])
            end
        end
        
        x_log = strcmp(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.isu_xu_parameteru_log.Text, 'true');
        y_log = strcmp(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{m}.isu_yu_parameteru_log.Text, 'true');
        
        % store information
        xmlGates.(gateName) = struct;
        xmlGates.(gateName).xparam = x_param;
        xmlGates.(gateName).yparam = y_param;
        xmlGates.(gateName).xlog = x_log;
        xmlGates.(gateName).ylog = y_log;
        xmlGates.(gateName).used = false;
    end
end

% check the consistency between xmlStruct and requirement (required_FACS_channels, clusterNum)
% 'P1' should be the FSC-SSC gate
if ~strcmp(xmlGates.P1.xparam, 'FSC-A') && ~strcmp(xmlGates.P1.yparam, 'SSC-A')
    error('P1 should be the FSC-SSC gate!')
end

% For clusterNum (e.g. 2) clusters to be sorted by numel(required_FACS_channels)
% channels (e.g. 3, channels A, B and C), 'P2' ~ 'P4' is AB, AC and BC,
% 'P5' ~ 'P7' is AB, AC and BC.
gateRecord = 2; % number of gate specified for this cluster
for n = 1:clusterNum
    used_channels = {};
    for m = 1:numel(required_FACS_channels)
        channelx = required_FACS_channels{m};
        for k = 1:numel(required_FACS_channels)
            channely = required_FACS_channels{k};
            if ~strcmp(channelx, channely) && isempty(find(contains(used_channels, channely), 1)) % combine two channels to have a gate
                thisGate = ['P', num2str(gateRecord)];
                if ~strcmp(xmlGates.(thisGate).xparam, channelx) || ~strcmp(xmlGates.(thisGate).yparam, channely)
                    error([thisGate, ' has incorrect channel setting!'])
                end
                gateRecord = gateRecord + 1;
            end
        end
        used_channels{end + 1} = channelx;
    end
end

end