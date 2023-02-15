function child=MultiMutationFunction1(positionindex,VarietyNumber,GenotypeLength,Pm,parent)
positionindex = 1:6;
if(rand<Pm && ~isempty(positionindex))
    position = [];
    for i = 1:length(positionindex)
        j = positionindex(i);
        position = [position (j-1)*sum(GenotypeLength)+1:(j-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3)+GenotypeLength(4)]; 
    end 
    index = randperm(length(position));
    mpoint = position(index(1:randi(length(position))));
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=[];
end