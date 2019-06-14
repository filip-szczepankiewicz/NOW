classdef optimizationProblem
    %OPTIMIZATIONPROBLEM Store all properties needed to run NOW.
    %   Parameters not specified by the user are default-initialized as follows:
    %
    %   Max gradient = 80 milliTesla/m
    %   Max slew rate = 100 milliTesla/m/milliSecond = 100 T/m/s
    %   Eta (heat dissipation parameter) = 1
    %   Discretization points = 77
    %   Target tensor = eye(3)
    %   Initialguess = 'random'
    %   zeroGradientAtIndex = [], i.e. only at start and end
    %   enforceSymmetry = false;
    %   redoIfFailed = true;
    %   useMaxNorm = false;
    %   doMaxwellComp = true;
    %   MaxwellIndex = 100;
    
    
    properties (Access = public)
        targetTensor = eye(3); % Isotropic encoding tensor
        N = 77;
        initialGuess = 'random';
        useMaxNorm = false;
        gMax = 80;
        sMax = 100;
        durationFirstPartRequested = 28;
        durationSecondPartRequested = 22;
        durationZeroGradientRequested = 8;
        eta = 1;
        enforceSymmetry = false;
        redoIfFailed = true;
        name = 'NOW';
        x0 = [];
        doMaxwellComp = true;
        doKmatrixNull = false;
        MaxwellIndex = 100;
        KmatrixIndex = 100;
        FlowIndex    = 50;
        AccIndex     = 50;
        MaxFunEval = 1e5;
        MaxIter    = 5e3;
    end
    
    properties (SetAccess = private)
        zeroGradientAtIndex = [];
        tolIsotropy = .5e-2; %before 1e-4
        tolMaxwell
        tolKmatrix
        tolFlow
        tolAcc
        signs
        tolSlew
        durationFirstPartActual
        durationZeroGradientActual
        durationSecondPartActual
        totalTimeActual
        dt
        gMaxConstraint
        sMaxConstraint
        integralConstraint
    end
    
    methods (Access = public)
        function obj = optimizationProblem(varargin)
            if nargin > 0
                settings = varargin{1};
                
                % Overwrite defaults with user-specified settings
                fieldNames = fieldnames(settings);
                for i = 1:length(fieldNames)
                    eval(['obj.' fieldNames{i} ' = getfield(settings, fieldNames{i});'])
                end
            end
            
            % Get actual times after discretization
            [obj.durationFirstPartActual, obj.durationZeroGradientActual, obj.durationSecondPartActual, obj.totalTimeActual, obj.zeroGradientAtIndex] = ...
                getActualTimings(obj.durationFirstPartRequested, obj.durationZeroGradientRequested, obj.durationSecondPartRequested, obj.N, obj.enforceSymmetry);
            
            
            % Compute private variables
            obj.dt = obj.totalTimeActual/obj.N; %Time step in milliseconds. Division by N instead of N-1 due to half step shift in gradients.
            obj.gMaxConstraint = obj.gMax*obj.dt;
            obj.sMaxConstraint = obj.sMax*obj.dt^2;
            obj.integralConstraint = obj.eta*obj.gMaxConstraint^2*obj.totalTimeActual/obj.dt;
            
            obj.tolMaxwell = obj.MaxwellIndex * obj.dt;
            obj.tolKmatrix = obj.KmatrixIndex * obj.dt;
            
            obj.tolFlow    = obj.FlowIndex / (obj.dt/1000)^4 /  (2.6751e+08)^2 * 1e6 * 1e6;
            obj.tolAcc     = obj.AccIndex  / (obj.dt/1000)^6 /  (2.6751e+08)^2 * 1e6 * 1e6;
            
            % Create sign vector to store info on spin dephasing direction.
            % This is now done independent of Maxwell compensation.
            if ~isempty(obj.zeroGradientAtIndex)
                signs = ones(obj.N - 1, 1); % Ghost points excluded during opt
                midpt = ceil(mean(obj.zeroGradientAtIndex));
                signs(midpt:end) = -1;
                obj.signs = signs;
            else
                obj.signs = ones(obj.N - 1, 1); % Ghost points excluded during opt
            end
            
            if obj.doMaxwellComp && isempty(obj.zeroGradientAtIndex)
                % Maxwell terms cannot be compensated if no 180 pulses are
                % present. Normally this kind of optimization is intended for
                % a repetition of two identical self-balanced waveforms,
                % but it may be necessary to warn users that single-sided
                % experiments will always incurr some error due to
                % concomitant fields. In practice the "weight" of the
                % optimization with respect to Maxwell terms is removed by
                % setting the tolerance to a large value.
                
                % This error can be commented out, or made into a warning, if you know 
                % what you are doing!
                error('Maxwell compensation is impossible for 0 refocusing pulses!')
            end
            
            if obj.doMaxwellComp == false
                obj.MaxwellIndex = 10^10;
            end
            
        end
        
    end
end


