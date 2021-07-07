function [] = plotPositions(Positions, posColors)

if ~isempty(Positions) && length(posColors) == length(Positions)+1
    hold on
    for j=1:length(Positions)
        drawpolygon('Position', Positions{j},'Color',posColors{j+1});
        hold on
        text(Positions{j}(3,1), Positions{j}(3,2), num2str(j),'FontSize',36,'Color','red')
        hold on
    end
end

xlabel('Phase V1');
ylabel('Phase V2');
title('REM')

end

