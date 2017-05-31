% EE239AS.2 Homework-3
% Utility function to plot spike trains
% Input: S is a vector of cells
%        each cell contains a vector of spike times 
%        eg. S{1} = [21, 39, 200, 659, 1999] 
%            S{2} = [1, 39, 387]

function plotRaster(S)

gap = 3;
mark = 5;
pad = 30;

numSpikeTrain = numel(S);

for s = 1:numSpikeTrain
    offset = pad + gap + (s-1)*(gap+mark);
    train = S{s};
    if ~isempty(train)
        for t = train'
            line([t t], [offset, offset+mark], 'Color', 'k');
        end
    end
end

xlabel('Time (ms)');
ylim([0 offset+mark+gap+pad]);
set(gca, 'YTickLabel', []);

end