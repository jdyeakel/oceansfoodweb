function nullweb(A)


    numlinks = sum(A);

    nullA = copy(A);

    #cartesian coordinates of links
    links = findall(!iszero,A);

    for i = 1:numlinks

        oldpos = links[i];

        #find empties
        emptylinks = findall(iszero,nullA);

        #choose new position
        newpos = emptylinks[rand(collect(1:length(emptylinks)))];

        #turn off old link; turn on new link
        nullA[oldpos] = 0;
        nullA[newpos] = 1;

    end

    return nullA
end