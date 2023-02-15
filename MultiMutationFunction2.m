function child=MultiMutationFunction2(VarietyNumber,GenotypeLength,Pm,parent)
if(rand<Pm)
    position = [];
    for i = 1:VarietyNumber
        position = [position (i-1)*sum(GenotypeLength)+1:(i-1)*sum(GenotypeLength)+GenotypeLength(1)+GenotypeLength(2)+GenotypeLength(3)+GenotypeLength(4)];
    end 
    index = randperm(length(position));
    mpoint = position(index(1:randi(length(position))));
    child = parent;
    child(mpoint) = 1-parent(mpoint);
else
    child=[];
end