function livingResult = livingSum(dead, nbirds, birdcur, axis)
    livingbirds = nbirds;
    value = 0;
    for i = 1:nbirds
        if birdcur(i,axis) == dead
            livingbirds = livingbirds - 1;
        else
            value = value + birdcur(i,axis);
        end
    end
    livingResult = value / livingbirds;
end