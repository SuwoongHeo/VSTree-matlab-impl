clear all; close all;

[pts, face] = icosphere(5);

vstree = VSTree(pts, 'ptsnormal', pts, 'nodeCapacity', 1, 'delta_d', 1/20);

colorid = unique(vstree.PointNodes);
colortab = repmat((colorid / max(colorid)),1,3).*rand(size(colorid,1), 3);
colortab = colortab./repmat(sqrt(sum(colortab.^2,2)),1,3);
color = zeros(size(vstree.PointNodes,1), 3);
for i=1:size(colorid,1)
    toDraw = vstree.PointNodes==colorid(i);
    color(toDraw,:) = repmat(colortab(i,:), sum(toDraw), 1);
end

figure, hd = vstree.plot(3, '--', 'Color', [0, 0.6, 1], 'linewidth', 0.8);
axis image, view(3);
hold on; dispMesh(pts, face, color, 1);