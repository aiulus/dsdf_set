% profile on
% <script_name>
% p = profile('info')
function customProfile(p, timeThreshold)
    % p: Profiling data structure from profile('info')
    % timeThreshold: Minimum execution time to include a function

    disp('Received profiling information.');

    % Step 1: Extract Call Hierarchy and Filter Based on Execution Time
    disp('Building filtered call hierarchy...');
    callList = {p.FunctionTable.FunctionName}; % All function names
    executionTimes = [p.FunctionTable.TotalTime]; % Execution times
    significantIndices = find(executionTimes >= timeThreshold); % Only keep significant functions

    % Create unique names for significant functions
    filteredCallList = arrayfun(@(idx) sprintf('%s [%d]', callList{idx}, idx), ...
                                significantIndices, 'UniformOutput', false);

    % Filter the edges for the graph
    edges = []; % Initialize edges for digraph
    for i = significantIndices
        childrenStruct = p.FunctionTable(i).Children; % Struct array for children
        if ~isempty(childrenStruct)
            childIndices = [childrenStruct.Index];
            % Keep only children that are also significant
            childIndices = intersect(childIndices, significantIndices);
            for j = 1:numel(childIndices)
                edges = [edges; find(significantIndices == i), find(significantIndices == childIndices(j))];
            end
        end
    end

    % Step 2: Create a Directed Graph
    disp('Visualizing filtered call hierarchy...');
    G = digraph(edges(:,1), edges(:,2), [], filteredCallList);

    % Step 3: Add Execution Times to Nodes
    filteredTimes = executionTimes(significantIndices); % Filter execution times
    nodeColors = filteredTimes / max(filteredTimes); % Normalize

    % Step 4: Plot the Graph
    figure;
    % h = plot(G, 'Layout', 'layered', 'NodeLabel', G.Nodes.Name, 'NodeCData', nodeColors);
    h = plot(G, 'Layout', 'layered', 'NodeLabel', G.Nodes.Name, 'NodeCData', nodeColors);
    h.XData = h.XData * 2; % Stretch horizontally
    h.YData = h.YData * 2; % Stretch vertically

    colorbar;
    title('Filtered Function Call Hierarchy');
    xlabel('Execution Time (Normalized)');
    set(h, 'EdgeColor', 'k'); % Set edge color for better visibility

    % Step 5: Save Results
    outputFile = sprintf('DDSF_filtered_call_hierarchy_threshold_%.3f.png', timeThreshold); % Dynamic file name
    disp(['Saving visualization to ', outputFile, '...']);
    saveas(gcf, outputFile); % Save the figure as an image
    disp('Filtered call hierarchy visualization saved.');
end
