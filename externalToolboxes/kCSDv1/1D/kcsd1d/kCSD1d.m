classdef kCSD1d < handle
    
    properties
        X         %estimation area
        src       %source positions
        
        R         %radius of disc for the kCSD1D method
        h         %width of basis element
        sigma
        
        lambda
        
        Rs
        
        lambdas
        
        csdEst         % estimated CSD
        
    end
    
    properties (SetAccess = private)
        elPos          % vector of electrode positions
        pots           % vector of potential values (n_electrodes, nt)
        distTable      % table used to interpolate potential value
        bPotMatrix
        KPot
        
        bSrcMatrix         % matrices for for estimating CSD & potentials
        interpCross
        bPot
        bSrc
        
        
    end
    
    
    %----------------------------------------------------------------------
    methods
        
        function k = kCSD1d(elPos, pots, varargin)
            k.elPos = elPos;
            k.pots = pots;
            
            
            propertyArgIn = varargin;
            while length(propertyArgIn) >= 2
                prop = propertyArgIn{1};
                val = propertyArgIn{2};
                propertyArgIn = propertyArgIn(3:end);
                switch prop
                    case 'X'
                        k.X = val;
                    case 'R'
                        k.R = val;
                    case 'h'
                        k.h = val;
                    case 'sigma'
                        k.sigma = val;
                    case 'gDx'
                        gDx = val;
                    case 'ext'
                        ext = val;
                    case 'xMin'
                        xMin = val;
                    case 'xMax'
                        xMax  = val;
                    case 'nSrc'
                        nSrc = val;
                        
                end
            end
            if ~exist('xMin', 'var') && ~exist('X', 'var')
                xMin = min(k.elPos(:,1));
            end
            if ~exist('xMax', 'var') && ~exist('X', 'var')
                xMax = min(k.elPos(:,1));
            end
            if ~exist('nSrc', 'var')
                nSrc = 300;
            end
            if ~exist('ext', 'var')
                ext = 0;
            end
            
            if isempty(k.X)
                k.X = meshgrid(xMin:gDx:xMax);
            end
            
            if isempty(k.R)
                k.R = 1;
            end
            
            if isempty(k.sigma)
                k.sigma = 1;
            end
            
            if isempty(k.h)
                k.h = elPos(2) - elPos(1);
            end;
            k.lambda = 0;
            k.src = makeSrc(k.X, ext, nSrc);
            dist_max = max(k.X(:)) - min(k.X(:));
            k.distTable = create_dist_table(100, dist_max, k.h, k.R, k.sigma, ...
                'gauss');
            k.bPotMatrix = bPotMatrixCalc(k.X, k.src, k.elPos, k.distTable);
            k.KPot=(k.bPotMatrix)'*(k.bPotMatrix);
            k.bSrcMatrix = bSrcMatrixCalc(k.X, k.src, k.h);
            k.interpCross = k.bSrcMatrix*k.bPotMatrix;
        end
        
        function estimate(k)
            [~, nt] = size(k.pots);
            nx = length(k.X);
            nEl = length(k.elPos);
            
            k.csdEst = zeros(nx, nt);
            KInv = (k.KPot + k.lambda.*eye(size(k.KPot)))^(-1);
            
            for t = 1:nt
                beta = KInv * k.pots(:,t);
                for i = 1:nEl
                    k.csdEst(:, t) = k.csdEst(:, t) + beta(i).*k.interpCross(:, i);
                end
            end
            
        end
        %----------------------------------------------------------------------
              
        
    end
    
end