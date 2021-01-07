function bsfoodweb(xmax,mass,trophic,occur)

    alpha = xmax[1];
    beta = xmax[2];
    gamma = xmax[3];

    
    ns = length(occur);
    nc = length(findall(x->x>1,trophic));
    amatrix = Array{Int64}(undef,ns,ns) .* 0;

    #Loop across consumers
    for i=1:ns
        if trophic[occur[i]] > 1
            # Loop across prey
            # The while loop means that all consumers will have at least one prey
            while sum(amatrix[:,i]) < 1
                for j=1:ns
                    #consumer mass
                    mi = mass[occur[i]];
                    #prey mass
                    mj =mass[occur[j]];

                    #The probability of an interaction
                    pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));

                    #draw
                    linkdraw = rand();
                    if linkdraw < pij
                        amatrix[j,i] = 1;
                    else
                        amatrix[j,i] = 0;
                    end
                end
            end
        end
    end

    #sort by body mass
    bm_sortperm = sortperm(mass[occur]);
    amatrix_sort = amatrix[bm_sortperm,bm_sortperm];

    return(amatrix_sort)

end

