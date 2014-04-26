function matFitness = lsdFitness(matrix)
xSize = size(matrix, 1);
ySize = size(matrix, 2);
[sampleX, sampleY] = find(matrix);
sampleCoords = [sampleX sampleY];
matFitness = sum(pdist(sampleCoords, 'euclidean'));
end