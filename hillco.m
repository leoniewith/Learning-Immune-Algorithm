function hillcoff = hillco(cellmem)
bandwidth=5;
[~,~,clustMembsCell] = MeanShiftCluster(cellmem,bandwidth);
celltemp = cellmem';
unicelltemp = celltemp;
hillcoff = [];
for i = 1:size(clustMembsCell,1)
    mcellindex = clustMembsCell{i};
    if length(mcellindex) == 1
        hillcoff(mcellindex) = 1;
    else
        clustertemp = unicelltemp(mcellindex,:);
        for j = 1:length(mcellindex)
            hillcoff(mcellindex(j)) = sum(max(clustertemp)-clustertemp(j,:)-(clustertemp(j,:)-min(clustertemp)));
        end
        if length(unique(hillcoff(mcellindex))) < length(hillcoff(mcellindex))
            [~,b] = max(hillcoff(mcellindex));
            hillcoff(mcellindex(b)) = hillcoff(mcellindex(b)) + 0.01;
        end
    end
end