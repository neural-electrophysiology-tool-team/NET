function mask=makeCircularMaskForGUI(radius,varargin)
% Create a single gaussian transparency mask and store it to a texture:
    interior_val=0;
    exterior_val=255;
    color = [1 0 0]*255;
%     pvpmod(varargin)

    sz=2*radius+1;
    mask=ones(sz, sz, 4);
    center = [floor(sz/2) floor(sz/2)];

    for r=1:sz
        for c=1:sz
            if sqrt((r-center(1))^2+(c-center(2))^2)<=radius
                mask(r,c,4)=interior_val;
            else
                mask(r,c,4)=exterior_val;
            end
        end
    end
    
    for i=1:3
        mask(:,:,i) = mask(:,:,i)*color(i);
    end   
end