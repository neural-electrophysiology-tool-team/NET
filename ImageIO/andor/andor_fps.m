function [fps, period, relerr] = andor_fps(namepattern, varargin)
%function [fps, period, relerr] = andor_fps(namepattern, varargin)
%calculates the fps and its inverse the inter-frame period of an andor sif
%file. Only the first frames are used for this calculation. They are
%assumed to be equidistant.
%
%input: filename (wildards allowed) of sif file
%optional input: maxframe: how many frames to use for the calculation (default 5)
%                maxerr: maximal relative error in equidistancy (default 0.02)
%
%output: fps, 
%        period, 
%        relerr, the relative error (jitter)



params.maxframe = 5;
params.maxerr = 0.02;
params = omex_read_params(params,varargin);

for ki = 1:params.maxframe
    try
        [data, meta] = sifread(wildfile(namepattern),ki);
        time(ki) = meta.TimeStampFrame / 1000000; %get time and convert from us to s
    catch ME %not all frames can be read
        if strcmpi(ME.identifier, 'check_nr_of_files:NoFile') || strcmpi(ME.identifier, 'check_nr_of_files:SeveralFiles');
            error(ME.message);
        else
            if ki ==1;
                error('andor_fps: Problem reading metadata. Is this really an Andor sif file?')
            else
                warning('andor_fps: Could not read as many frames as requested for mean. I will do with %d frames', ki-1);
            end
        end
        break
    end
end

dt = diff(time);

if isempty(dt)%we have a singel image and no difference was calculated.
    period = Inf;
    fps = 0;
    relerr = NaN;
    warning('andor_fps: Only one frame. No framerate can be calculated');
else
    [period, fps, dt, relerr ] = goodmean(dt);
end

if isnan(period)% if the error is large, try again without hte first frame (which was deleted by the goodmean function
    [period, fps, dt, relerr ] = goodmean(dt);
end



    %subfunction to calculate mean and check if it is reliable:
    
    function [period, fps, dt, relerr] = goodmean(dt)
        dtmean = mean(dt);
        dtstd = std(dt);
        relerr = dtstd/dtmean;
        
        if relerr < params.maxerr
            period = dtmean;
            fps = 1/period;
        else %the first frame might be different, remove it, then we can try again 
            period = NaN;
            fps = NaN;
            dt(1) = [];
        end
    end

end
