function editxml2(xmlCluster, xmlFile, channelsOfInterest, output_dir)

% This function is to edit the xmlFile template from FACS machine,
% according to the xmlCluster calculated by CART tree and gateID. The
% xmlFile is firstly read as a struct, edited and finally stoed as a .xml
% file that is readily available to be used in cell sorting.
% 
% Usage: editxml(xmlCluster, xmlFile, channelsOfInterest, output_dir)
%
% Input
% The input xmlCluster is a struct mimicking the xml file, describing two 4-point
% polygon gates used to obtain each leaf node belonging to one cluster. In particular, it is a
% hierarchical struct with structures like the following.
% ClusterX -> a cell array, in this cell array, we have a hierarchical
% struct, with two fields: gates and properties.
% Firstly, gates{gate1, gate2, ...} (gate is a struct) -> xparam, yparam,
% used (whether it will be used to write in the final .xml file), 
% points{point1, point2, point3, point4} (point is a struct) -> x, y
% Thus, you can access the first point of gate1 of cluster 1 using:
% xmlCluster.Cluster1.gates{1}.points{1}.x(or y)
% Note that here is just the gates calculated from CART tree. It can be
% further refined by additional processing.
% Secondly, properties is a struct, pointing to three fields, purity,
% relativeYield and absoluteYield. Purity simply represents the purity for
% this leaf node; relativeYield represents the yield of this leaf node
% relative to the cells in the relevant cluster; absoluteYield represents
% the yield of this leaf node relative to all the cells.
% See saveTree2 for details.
% 
% The input xmlFile is a .xml file from FACS machine, with standard
% Gating-ML format. The xmlFile should have nchoosek(channelNum,
% 2)*clusterNum + 1 (x:FSC, y:SSC) + 1 (All Events) gates (e.g. 4 channels, 5 clusters, thus
% 6*5 + 1 + 1 = 32 gates). Some of them will be editted according to xmlCluster. Note
% that we expect that each gate is a 4-points-polygon defined in 2D space,
% with x and y axis specifying one channel, respectively. In particular,
% fluorescent channel of x axis should always has excitation wavelength
% smaller than that of y. For instance, x = 'Pacific Blue-A' and y =
% 'FITC-A' are legal and vice versa is not. Note that 'All Events' is
% parent of FSC-SSC gate; and FSC-SSC gate is parent to all fluorescent gate.
% In particular, given three fluorescent channels (A, B, C) and 2 cell subpopulations\
% to be sorted, 'P1' should be FSC-SSC gate, 'P2' ~ 'P4' is AB, AC and BC,
% 'P5' ~ 'P7' is AB, AC and BC.
%
% The input channelsOfInterest is a cell array, specifying the channels of
% interest. Accepted channels include 'CFP', 'GFP', 'mCherry' and 'APC'
%
% The input output_dir is the integrate path for storage of the output
% files, such as 'E:\Data\D0500\'.
%
% Output
% According to the path specified by output_dir, two files will be saved.
% The first one is the edited final_gate.xml file, used in subsequent cell
% sorting; the second is multiple gate2cluster.csv files, storing how
% combinations of gates lead to each desired leaf node belonging to one
% particular cluster. Note that for each leaf node, one or two gates (4-points-polygon in 2D space) are combined (&&) to
% obtain the desired cell subpopulation. Each file will be named after as
% below: "ClusterX_gate2cluster.csv"
%
% Dependency: xmlGateExtract.m by Tianmin Wang, Tsinghua
%             xml2struct by W. Falkena, ASTI, TUDelft
% 
% Written by Ping Shen and Tianmin Wang @ Tsinghua University
% Version 0.2.2. Created on Oct 19, 2019. Last modified on May 14, 2021.

% struct storing common facs names and relation to microscopy channel name
facs_mapping = struct;
facs_mapping.FSC = 'FSC-A';
facs_mapping.SSC = 'SSC-A';
facs_mapping.CFP = 'Pacific Blue-A';
facs_mapping.GFP = 'FITC-A';
facs_mapping.mCherry = 'PE-Texas Red-A';
facs_mapping.APC = 'APC-A';

% convert channels to FACS channels
required_FACS_channels = {};
for m = 1:numel(channelsOfInterest)
    channel = channelsOfInterest{m};
    required_FACS_channels{end + 1} = facs_mapping.(channel);
end

% check the consistency between xmlCluster and xmlFile
clusterNum = numel(fields(xmlCluster));
xmlStruct = xml2struct(xmlFile); % read in .xml file as a struct
xmlGates = xmlGateExtract(xmlStruct, required_FACS_channels, clusterNum); % extract information about all gates, and check consistency
allGates = fields(xmlGates);

% initialization of cluster_gate_association, used to save the gating
% strategy for each cluster
all_clusters = {};
cluster_gate_association = struct; % ClusterX -> a table with n rows (leaf node) and 3 columns (leaf node and two gates)
for m = 1:clusterNum
    cluster = ['Cluster', num2str(m)];
    all_clusters{end + 1} = cluster;
    n = numel(xmlCluster.(cluster)); % number of leaf node for this cluster
    cluster_gate_association.(cluster) = cell2table((num2cell(1:n))', 'VariableNames',{'LeafNode'});
    gate_initialization = {};
    property_initialization = zeros(numel(xmlCluster.(cluster)), 1);
    for k = 1:numel(xmlCluster.(cluster))
        gate_initialization{end + 1} = '';
    end
    cluster_gate_association.(cluster).gate1 = gate_initialization';
    cluster_gate_association.(cluster).gate2 = gate_initialization';
    cluster_gate_association.(cluster).purity = property_initialization;
    cluster_gate_association.(cluster).relativeYield = property_initialization;
    cluster_gate_association.(cluster).absoluteYield = property_initialization;
end

% edit the xmlStruct to incorporate the xmlCluster about calculated gate
% coordinations
for m = 1:clusterNum
    cluster = all_clusters{m};
    for h = 1:numel(xmlCluster.(cluster))
        % write leaf node purity and yield into cluster_gate_association
        cluster_gate_association.(cluster)(h, 'purity') = num2cell(xmlCluster.(cluster){h}.properties.purity);
        cluster_gate_association.(cluster)(h, 'relativeYield') = num2cell(xmlCluster.(cluster){h}.properties.relativeYield);
        cluster_gate_association.(cluster)(h, 'absoluteYield') = num2cell(xmlCluster.(cluster){h}.properties.absoluteYield);
        for k = 1:numel(xmlCluster.(cluster){h}.gates)
            if xmlCluster.(cluster){h}.gates{k}.active % active gate
                % get the metric about this gate
                xparam = xmlCluster.(cluster){h}.gates{k}.xparam;
                xparam = facs_mapping.(xparam);
                yparam = xmlCluster.(cluster){h}.gates{k}.yparam;
                yparam = facs_mapping.(yparam);
                % bypassing negative values and consequent problem once log10 rescaling is needed
                point1_x = max(xmlCluster.(cluster){h}.gates{k}.points{1}.x, 1);
                point1_y = max(xmlCluster.(cluster){h}.gates{k}.points{1}.y, 1);
                point2_x = max(xmlCluster.(cluster){h}.gates{k}.points{2}.x, 1);
                point2_y = max(xmlCluster.(cluster){h}.gates{k}.points{2}.y, 1);
                point3_x = max(xmlCluster.(cluster){h}.gates{k}.points{3}.x, 1);
                point3_y = max(xmlCluster.(cluster){h}.gates{k}.points{3}.y, 1);
                point4_x = max(xmlCluster.(cluster){h}.gates{k}.points{4}.x, 1);
                point4_y = max(xmlCluster.(cluster){h}.gates{k}.points{4}.y, 1);
            
                % edit xml
                gateAvailabilityFlag = false;
                for n = 1:numel(allGates)
                    gate = allGates{n};
                    gate_xchannel = xmlGates.(gate).xparam;
                    gate_ychannel = xmlGates.(gate).yparam;
                    if strcmp(gate_xchannel, xparam) && strcmp(gate_ychannel, yparam) && ~xmlGates.(gate).used % correct gate
                        
                        % rescale if necessary
                        if xmlGates.(gate).xlog
                            point1_x = log10(point1_x);
                            point2_x = log10(point2_x);
                            point3_x = log10(point3_x);
                            point4_x = log10(point4_x);
                        end
                    
                        if xmlGates.(gate).ylog
                            point1_y = log10(point1_y);
                            point2_y = log10(point2_y);
                            point3_y = log10(point3_y);
                            point4_y = log10(point4_y);
                        end
                    
                        % find the relevant gate in xmlStruct to write
                        for p = 1:numel(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate)
                            if strcmp(xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.name.Text, gate) % find the relevant gate in xml struct
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{1}.Attributes.x = num2str(point1_x);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{1}.Attributes.y = num2str(point1_y);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{2}.Attributes.x = num2str(point2_x);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{2}.Attributes.y = num2str(point2_y);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{3}.Attributes.x = num2str(point3_x);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{3}.Attributes.y = num2str(point3_y);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{4}.Attributes.x = num2str(point4_x);
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.region.points.point{4}.Attributes.y = num2str(point4_y);
                                xmlGates.(gate).used = true;
                                gateAvailabilityFlag = true;
                                break;
                            end
                        end
                    
                        % assign the relevant to its cluster, saved in cluster_gate_association
                        if xmlGates.(gate).used
                            gate1_state = table2cell(cluster_gate_association.(cluster)(h, 'gate1'));
                            gate1_state = gate1_state{1};
                            gate2_state = table2cell(cluster_gate_association.(cluster)(h, 'gate2'));
                            gate2_state = gate2_state{1};
                            if strcmp(gate1_state, '') % gate1 of this cluster is empty
                                cluster_gate_association.(cluster)(h, 'gate1') = {gate};
                            elseif strcmp(gate2_state, '') % gate1 has been assigned; gate2 is available
                                cluster_gate_association.(cluster)(h, 'gate2') = {gate};
                                original_parent = xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.parent.Text;
                                xmlStruct.bdfacs.experiment.acquisitionu_worksheets.worksheetu_template.gates.gate{p}.parent.Text = [original_parent, '\', gate1_state];
                            else
                                error(['more than two gates are active for ', cluster, '!'])
                            end
                            break;
                        end
                        
                    end
                end
                if ~gateAvailabilityFlag
                    error(['The number of ', xparam, ' vs. ', yparam, ' gates in template .xml file is not enough!\n\n'])
                end
            end
        end
    end
end

% generate the .xml file
mystruct2xml(xmlStruct, [output_dir, '\calculated_gate.xml']);

% save the result about gate combinations that lead to different clusters
% N .csv files will be saved, corresponding to each cluster.
for m = 1:clusterNum
    cluster = ['Cluster', num2str(m)];
    writetable(cluster_gate_association.(cluster), [cluster, '_final_gating_strategy.csv'], 'Delimiter',',');
end