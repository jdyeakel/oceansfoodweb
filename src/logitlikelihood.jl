function logitlikelihood(A,massvec,alphavec,betavec,gammavec)


    # alphavec = collect(-10:1:10);
    # betavec = collect(-10:1:10);
    # gammavec = collect(-10:1:10);
    numprey = size(A)[1];
    numpred = size(A)[2];
    
    massvec_pred = massvec;
    massvec_prey = massvec;
    
    L = SharedArray{Float64}(length(alphavec),length(betavec),length(gammavec));
    nparam = length(alphavec)*length(betavec)*length(gammavec);
    avec = repeat(collect(1:length(alphavec)),outer=length(betavec)*length(gammavec),inner=1);
    bvec = repeat(collect(1:length(betavec)),inner=length(alphavec),outer=length(gammavec));
    gvec = repeat(collect(1:length(gammavec)),outer=1,inner=length(alphavec)*length(betavec));
    
    # avec = repeat(collect(1:length(alphavec)),inner=length(alphavec))
    
    
    stepmatrix = [avec bvec gvec];

    @sync @distributed for k = 1:nparam
        a = stepmatrix[k,1];
        b = stepmatrix[k,2];
        g = stepmatrix[k,3];
        
        alpha = alphavec[a];
        beta = betavec[b];
        gamma = gammavec[g];
        
        let cumLij = 0.0
            for i=1:numpred
                for j=1:numprey
                    if i != j
                        aij = copy(A[j,i]);
                        mi = massvec_pred[i]; #i is predator
                        mj = massvec_prey[j]; #j is prey
                        
                        Lij = 0.0;
                        
                        pij = exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2)/(1+exp(alpha + beta*log(mi/mj) + gamma*log(mi/mj)^2));
                        
                        if aij == 1 #ERROR: 7/15/20 THIS WAS A ZERO HA
                            Lij = copy(pij);
                        else
                            Lij = copy(1 - pij);
                        end
                        
                        cumLij += log(Lij);
                    end
                end
            end
            L[a,b,g] = copy(cumLij);
        end
        
    end
    
    maxL = findmax(L);
    thetavec = [alphavec[maxL[2][1]],betavec[maxL[2][2]],gammavec[maxL[2][3]]];
    
    alphamax = thetavec[1];
    betamax = thetavec[2];
    gammamax = thetavec[3];
    
    # global cumnum = 0.0;
    # global cumdenom = 0.0;
    cnumerator = Array{Float64}(undef,1);
    cdenominator = Array{Float64}(undef,1);
    let cumnum = 0.0, cumdenom = 0.0
        for i=1:numpred
            for j=1:numprey

                aij = copy(A[j,i]);
                mi = massvec_pred[i]; #i is predator
                mj = massvec_prey[j]; #j is prey
                
                pij = exp(alphamax + betamax*log(mi/mj) + gammamax*log(mi/mj)^2)/(1+exp(alphamax + betamax*log(mi/mj) + gammamax*log(mi/mj)^2));
                
                num = aij*pij + (1-aij)*(1-pij);
                denom = aij + (1-aij);
                
                
                cumnum += num;
                cumdenom += denom;
            end
        end
        cnumerator[1] = copy(cumnum);
        cdenominator[1] = copy(cumdenom);
    end
    
    fc = cnumerator[1]/cdenominator[1];
    
    return L, thetavec, fc
    
end

                