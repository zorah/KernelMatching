function [cells, indices] = fast_voronoi(X, X_lores, n_cells, seeds, samples, mesh)

if n_cells == numel(seeds)
    centroids = seeds;
else
    centroids = fps_euclidean(X.VERT(seeds,:), n_cells, 1);
    centroids = seeds(centroids);
end

remaining_samples = [setdiff(seeds, centroids); samples];
[centroids_lores,~] = knnsearch(X_lores.VERT, X.VERT(centroids,:));
[remaining_samples_lores,~] = knnsearch(X_lores.VERT, X.VERT(remaining_samples,:));

inside = false;
if ~isfield(X_lores, 'D')
    if ~exist('mesh', 'var')
        inside = true;
        try
            mesh = geodesic_new_mesh(X_lores.VERT, X_lores.TRIV);
        catch ex
            geodesic_delete;
            rethrow(ex);
        end
    end
    
    algorithm = geodesic_new_algorithm(mesh, 'dijkstra');
    
    source_points = cell(1,length(centroids));
    for i=1:length(source_points)
        source_points{i} = geodesic_create_surface_point(...
            'vertex',centroids_lores(i),X_lores.VERT(centroids_lores(i),:));
    end

    geodesic_propagate(algorithm, source_points);


cells = zeros(length(remaining_samples),1);
for i=1:length(remaining_samples)
        q = geodesic_create_surface_point('vertex', remaining_samples_lores(i), X_lores.VERT(remaining_samples_lores(i),:));
        [cells(i), ~] = geodesic_distance_and_source(algorithm, q);
end
if inside
    geodesic_delete;
end

else
    [~, cells] = min(X_lores.D(centroids_lores,remaining_samples_lores));
end

%     samples = [centroids; remaining_samples];
cells = [(1:numel(centroids))' ; cells];
indices = [centroids; remaining_samples];
end

