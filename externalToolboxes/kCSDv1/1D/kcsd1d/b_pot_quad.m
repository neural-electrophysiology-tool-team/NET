function [y]=b_pot_quad(src,arg,h,R,sigma,src_type)
switch src_type
    case 'step'
        y=quad ( @(zp) pot_intarg(zp,arg,h,R,sigma,src_type), src-(0.5*h),src+(0.5*h));
    case 'gauss'
        y=quad ( @(current_pos) pot_intarg(src, arg, current_pos, h, R, sigma, src_type), src-4*h,src+4*h);
end;


