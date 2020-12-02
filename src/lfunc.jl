function lfunc(x,A,massvec)
    #Calculate Likelihood of proposed alpha, beta, gamma
    alpha = x[1];
    beta = x[2];
    gamma = x[3];
    
    numprey = size(A)[1];
    numpred = size(A)[2];
    massvec_pred = massvec;
    massvec_prey = massvec;
    
    L = Array{Float64}(undef,1);
    
    let cumLij = 0.0
        for i=1:numpred
            for j=1:numprey
                # if i != j
                    
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
                # end
            end
        end
        L[1] = -cumLij
    end
    return L[1]
end