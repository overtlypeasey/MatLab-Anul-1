classdef DataManager
    properties
        DataFolderPath = fullfile(pwd); % Default data
        CurrentData % Stores the currently loaded data
    end

    methods 
        function fileList = listDataFiles(obj)
            % List all data files in the data folder
            d = dir(fullfile(obj.DataFolderPath, '*.mat'));
            fileList = {d.name};
        end
        
        function obj = loadData(obj, fileName)
            % Load data from a file
            loadedData = load(fullfile(obj.DataFolderPath, fileName));
            obj.CurrentData = loadedData;
            
        end
        function processedData = processData(obj, plotType)
            % Process data based on the selected plot type
            % This is a placeholder for actual data processing
            processedData = obj.CurrentData; % Placeholder
        end
    end
end