function saveLegendToImage(figHandle, legHandle, fileName, fileType)

%make all contents in figure invisible
allLineHandles = findall(figHandle, 'type', 'line');

for i = 1:length(allLineHandles)

    allLineHandles(i).XData = NaN; %ignore warnings

end

%make axes invisible
axis off

%move legend to lower left corner of figure window
legHandle.Units = 'pixels';
boxLineWidth = legHandle.LineWidth;
%save isn't accurate and would swallow part of the box without factors
legHandle.Position = [6 * boxLineWidth, 6 * boxLineWidth, ...
    legHandle.Position(3), legHandle.Position(4)];
legLocPixels = legHandle.Position;

%make figure window fit legend
%figHandle.Units = 'pixels';
figHandle.InnerPosition = [1, 1, legLocPixels(3) + 12 * boxLineWidth, ...
    legLocPixels(4) + 12 * boxLineWidth];

%save legend
saveas(figHandle, [fileName, '.', fileType], fileType);

end