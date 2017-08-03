function [h] = vis_partial_matches( X, Y, xin, yin)
%VIS_PARTIAL_MATCHES Summary of this function goes here
%   Detailed explanation goes here

colorsY = create_colormap(Y,Y);
colorsY = colorsY(yin,:);
h = figure;
subplot(1,2,1),
trisurf(Y.TRIV, Y.VERT(:,1), Y.VERT(:,2), Y.VERT(:,3), ones(Y.n, 1), 'EdgeAlpha', 0), hold on,
lighting phong, camlight,
scatter3(Y.VERT(:,1), Y.VERT(:,2), Y.VERT(:,3), 30, 0.8*ones(Y.n, 3), 'filled'),
axis equal, axis off,
hold on,
scatter3(Y.VERT(yin,1), Y.VERT(yin,2), Y.VERT(yin,3), 80, colorsY, 'filled'),
subplot(1,2,2),
white = 200 * [1 1 1];
colormap white,
trisurf(X.TRIV, X.VERT(:,1), X.VERT(:,2), X.VERT(:,3), ones(X.n, 1), 'EdgeAlpha', 0), hold on,
lighting phong, camlight,
axis equal, axis off,
scatter3(X.VERT(:,1), X.VERT(:,2), X.VERT(:,3), 50, 0.8*ones(1,3), 'filled'), hold on
scatter3(X.VERT(xin,1), X.VERT(xin,2), X.VERT(xin,3), 80, colorsY, 'filled'), 

end

