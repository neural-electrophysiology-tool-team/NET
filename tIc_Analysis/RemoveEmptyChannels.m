function [newT,newic]=RemoveEmptyChannels(t,ic)
%function [newT,newic]=RemoveEmptyChannels(t,ic);
% a function that removes empty channels data 

x=[];
for i=1:size(ic,2)
    if ~(ic(4,i)>=ic(3,i))
        x=[x ic([3 4],i)];
    end
end
if isempty(x)
    newT=t;
    newic=ic;
    disp('No empty channels were detected');
    return;
end

for i=1:size(x,2),
    fprintf('removing neuron [%d;%d]\n',x(1,i),x(2,i));
    index=find(ic(1,:)==x(1,i) & ic(2,:)==x(2,i));
    %cut out a peace of t:
    cutT=t(ic(3,index):ic(4,index));
    if ic(3,index)~=1,
        tBefore=t(1:ic(3,index)-1);
    else
        tBefore=[];
    end
    if ic(4,index)~=length(t),
        tAfter=t(ic(4,index)+1:end);
    else
        tAfter=[];
    end
    t=[tBefore,tAfter];
    %update indexes:
    ic(3,find(ic(3,:)>ic(3,index)))=ic(3,find(ic(3,:)>ic(3,index)))-length(cutT);
    ic(4,find(ic(4,:)>ic(4,index)))=ic(4,find(ic(4,:)>ic(4,index)))-length(cutT);
end
%remove channels:
for i=1:size(x,2),
        indexes(i)=find(ic(1,:)==x(1,i) & ic(2,:)==x(2,i));
end
ic(:,indexes)=[];

newT=t;
newic=ic;

%%%%%%% To remove neurons that dont fire (instead of the nan assigned to them);
for i=1:size(newic,2)
    if isnan(newT(newic(3,i)))
        newT(newic(3,i))=[];
        newic(:,i)=NaN;
        newic(3:4,(i+1):end)=newic(3:4,(i+1):end)-1;
    end
end
newic(:,isnan(newic(1,:)))=[];        
    

