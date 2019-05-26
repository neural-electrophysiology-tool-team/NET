function mask=makeRectangularMaskForGUI(txtSrectWidth,txtSrectHeight,varargin)
% Create a triangular transparency mask and store it to a texture:

try
    interior_val=0;
    exterior_val=255;
    gray=0.5;
%     pvpmod(varargin)
    
    x = txtSrectWidth/2;
    y = txtSrectHeight/2;
    
    sz = max(txtSrectWidth,txtSrectHeight);
    mask = ones(sz, sz, 2) * gray;
    center = [floor(sz/2) floor(sz/2)];
    
    for r=1:sz
        for c=1:sz
            if abs(r-center(1))<=x && abs(r-center(2))<=y
                mask(r,c,2)=interior_val;
            else
                mask(r,c,2)=exterior_val;
            end
        end
    end
catch ME
    rethrow(ME);
end