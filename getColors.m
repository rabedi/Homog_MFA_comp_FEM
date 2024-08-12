function lc = getColors(colorNum, isModel, includeBlack)

if (isModel)
    lc{1} = 'b';
    lc{2} = 'r';
    lc{3} = 'k';
    lc{4} = 'g';
    return;
end
lc = cell(0);
cntr = 1;
if (includeBlack)
    lc{cntr} = [0 0 0]; % black
    cntr = cntr + 1;
end
lc{cntr} = [1	0	0]; %red'
cntr = cntr + 1;
lc{cntr} = [1	102/255	0]; %orange'
cntr = cntr + 1;
lc{cntr} = [0	0	1]; % blue 'b';
cntr = cntr + 1;
lc{cntr} = [0 135/255 0]; %green
cntr = cntr + 1;
lc{cntr} = [0	1	1]; % teal [0	0	0.5]; % dark blue
cntr = cntr + 1;
lc{cntr} = [0.5	0.25	0]; % brown
%lc = {'b', 'r', 'k', 'g'};
cntr = cntr + 1;
lc{cntr} = 1/255*[255	0	255]; % magenta
cntr = cntr + 1;
lc{cntr} = 1/255*[255	128	192]; % 'rosy_pink'
cntr = cntr + 1;
lc{cntr} = [0.5	0.5	0.5]; %'dark_gray2'
cntr = cntr + 1;
lc{cntr} = [0	0.5	0.25]; %green blue
cntr = cntr + 1;
lc{cntr} = [0.75	0.75	0.75]; %'gray2'
cntr = cntr + 1;
lc{cntr} = [1	1	0]; %yellow , [0	0.5	0.25]; %green blue
cntr = cntr + 1;
lc{cntr} = [0.5	0	1]; % purple 
cntr = cntr + 1;
lc{cntr} = 1/255 * [203 0   51]; % red2
cntr = cntr + 1;
lc{cntr} = [0.5	0.5	0]; % olive 
cntr = cntr + 1;
lc{cntr} = 1/255 * [64	128	128]; % blue2
cntr = cntr + 1;
lc{cntr} = 1/255 * [255	128	192]; % rosy_pink
cntr = cntr + 1;
lc{cntr} = 1/255 * [255	128	128]; % peach
cntr = cntr + 1;
lc{cntr} = 1/255 * [128	0	64]; % arghavani
cntr = cntr + 1;
lc{cntr} = 1/255 * [128	0	128]; % purple2
szc = length(lc);
for i = szc + 1:99
    ii = mod(i - 1, szc) + 1;
    lc{i} = lc{ii};
end

