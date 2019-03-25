rc=atsif_setfileaccessmode(0);
rc=atsif_readfromfile('spectrum.sif');
if (rc == 22002)
  signal=0;
  [rc,present]=atsif_isdatasourcepresent(signal);
  if present
    [rc,no_frames]=atsif_getnumberframes(signal);
    if (no_frames > 0)
        [rc,size]=atsif_getframesize(signal);
        [rc,left,bottom,right,top,hBin,vBin]=atsif_getsubimageinfo(signal,0);
        xaxis=0;
        [rc,data]=atsif_getframe(signal,0,size);
        [rc,pattern]=atsif_getpropertyvalue(signal,'ReadPattern');
        if(pattern == '0')
           calibvals = zeros(1,size);
           for i=1:size,[rc,calibvals(i)]=atsif_getpixelcalibration(signal,xaxis,(i)); 
           end 
           plot(calibvals,data);      
           title('spectrum');
           [rc,xtype]=atsif_getpropertyvalue(signal,'XAxisType');
           [rc,xunit]=atsif_getpropertyvalue(signal,'XAxisUnit');
           [rc,ytype]=atsif_getpropertyvalue(signal,'YAxisType');
           [rc,yunit]=atsif_getpropertyvalue(signal,'YAxisUnit');
           xlabel({xtype;xunit});
           ylabel({ytype;yunit});
        elseif(pattern == '4')
            width = ((right - left)+1)/hBin;
            height = ((top-bottom)+1)/vBin;
           newdata=reshape(data,width,height);
           imagesc(newdata);
        end
    end    
  end
  atsif_closefile;
else
  disp('Could not load file.  ERROR - ');
  disp(rc);
end
