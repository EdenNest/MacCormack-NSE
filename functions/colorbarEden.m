function cb=colorbarEden(str,varargin)
    % This is my own function to clean up visualization. Colorbar with
    % latex is messy :) It doesn't define Latex as main interpretor (as of
    % Feb 2024) with the general groot command

    % varargin is to set limits (arbitrary). str should be in Latex.

    cb=colorbar;
    cb.TickLabelInterpreter = 'latex';
    cb.Label.Interpreter = 'latex';
    cb.Label.String = str;
    switch nargin
        case 1
        case 3
            clim([varargin{1} varargin{2}]);
        otherwise
            disp("wrong colorbarEden input")
    end
end