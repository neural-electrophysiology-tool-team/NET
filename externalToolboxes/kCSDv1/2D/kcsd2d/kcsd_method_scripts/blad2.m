function er = blad2(data1, data2)
    alpha = sum(squeeze(data1(:)).*data2(:))...
        /sum(squeeze(data1(:)).^2);
    er = sum((data1(:)*alpha-data2(:)).^2)/sum(data2(:).^2);