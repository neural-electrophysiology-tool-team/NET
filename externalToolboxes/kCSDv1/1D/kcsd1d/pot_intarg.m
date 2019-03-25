function [y]=pot_intarg(src, arg, current_pos, h, R, sigma, src_type)

        y=(1/(2*sigma))*(sqrt((arg-current_pos).^2+R.^2)-abs(arg-current_pos)) .* ...
            gauss_rescale(src, current_pos, h);     
end