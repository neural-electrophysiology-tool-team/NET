% Function purpose : Returns vectors that represent binary words of
% activation (groups of neurons firing together)
%
% Function recives :    t [ms] - firing times
%                       ic - index channel
%                       triggers - event times
%                       win - the length of each event
%                       res - the resolution of the SB detection (bin size)
%                       sigma - gaussian width in the convolution
%                       mediaWindow - the window for the floating median
%
% Function give back :  binaryVectors - [numberOfEvents x numberOfNeurons]
%                       contains binary words of activation, for each small 
%                       event the eventLocator identified.
%                       intensities [input voltage units] - the intensities
%                       for each activation
%                       durations [ms] - the durations of each activation
%                       (the resolution is defined by res)
% 
function [binaryVectors, intensities, durations]=getBinaryVectors(t, ic, triggers, win, res, sigma, medianWindow)
    binaryVectors = [];
    intensities = [];
    durations = [];
    eventNumber = 0;
    h = waitbar(0,'Extracting binary vectors');
    for eventStart = triggers'
        eventNumber = eventNumber + 1;
        [tN,icN]=CutSortChannel(t,ic,eventStart,eventStart+win,0);

        [BS,~,BE,BI,~,~,~,~,BD]=eventLocator(ones(1,length(tN)),tN-eventStart,icN,20,'res',res,'sigma',sigma,'medianWindow',medianWindow,'addToSides',15,'stdThresh',1.5,'plotResults',0);
        if isempty(BS) || isempty(BE)
            fprintf("Event number %d was empty\n", eventNumber);
            continue
        end
        M=BuildBurstMatrix(ic,round(t/res),round(eventStart/res),round(win/res));
        M(M>1)=1;
        
        m = squeeze(M);
        eventTimes = [floor(BS/res);floor(BE/res)];
        mBinaryVectors = zeros([size(eventTimes,2),361]);
        i = 0;
        for event = eventTimes
            i=i+1;
            mBinaryVectors(i,:) = any(m(:,event(1):event(2))');
        end
        
        binaryVectors = [binaryVectors; mBinaryVectors];
        intensities = [intensities; BI'];
        durations = [durations; BD'];
        
        waitbar(eventNumber / size(triggers',2))
    end
    close(h)
end