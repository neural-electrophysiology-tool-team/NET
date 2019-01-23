load('layout_40_16x2_FlexLin.mat');

padSize=[12 35];
hight=600;
width=105;
dx=width/15;
dy=hight/15;

linearOrder=[24 25 27 26 22 23 28 29 20 21 31 30 18 19 1 32 16 17 14 15 3 2 12 13 4 5 10 11 7 6 8 9];

tmpPos=[0;0];
for i=1:16
    pElec=find(linearOrder==i);
    Enp(:,pElec)=tmpPos+[dx;-dy];
    tmpPos=Enp(:,pElec);
end

tmpPos(1)=tmpPos(1)+15-dx;
tmpPos(2)=tmpPos(2)-dy;
for i=17:32
    pElec=find(linearOrder==i);
    Enp(:,pElec)=tmpPos+[dx;dy];
    tmpPos=Enp(:,pElec);
end
Enp=Enp-min(Enp')'*ones(1,32);

save('layout_40_16x2_FlexLin.mat','En','Ena','Enp');


