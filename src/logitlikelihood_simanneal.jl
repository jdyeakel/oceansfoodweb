function logitlikelihood(A,massvec,initial)
    
    L = Array{Float64}(undef,0);
    push!(L,-10000);
    
    alphavec = Array{Float64}(undef,0);
    betavec = Array{Float64}(undef,0);
    gammavec = Array{Float64}(undef,0);
    tvec = Array{Float64}(undef,0);
    
    let deltaL = 1.0, 
        newalpha = initial[1],
        newbeta = initial[2],
        newgamma = initial[3],
        tic = 0,
        temp = 100.
        
        deltaL = 1.0; newalpha = initial[1]; newbeta = initial[2]; newgamma = initial[3]; tic = 0; temp = 100.
        
        while tic < 10000
            
            #Starting values of alpha, beta, gamma
            oldalpha = copy(newalpha);
            oldbeta = copy(newbeta);
            oldgamma = copy(newgamma);
            
            #Propose an altered alpha, beta, gamma
            if tic > 0
                alpha = newalpha + newalpha*rand(Normal(0,0.01));
                beta = newbeta + newbeta*rand(Normal(0,0.01));
                gamma = newgamma + newgamma*rand(Normal(0,0.01));
            else
                alpha = copy(newalpha);
                beta = copy(newbeta);
                gamma = copy(newgamma);
            end
            
            #Calculate Likelihood of proposed alpha, beta, gamma
            numprey = size(A)[1];
            numpred = size(A)[2];
            massvec_pred = massvec;
            massvec_prey = massvec;
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
                global L = push!(L,cumLij);
            end
            
            #Do we accept this new alpha, beta, gamma?
            #Calculate cceptance probability
            paccept = 0.0;
            if L[end] > L[end-1]
                paccept = 1.0;
            else
                if L[end] != -Inf
                    paccept = exp(-(L[end-1] - L[end])/temp);
                else 
                    paccept = 0.0;
                end
                # paccept = 0.0;
            end
            pdraw = rand();
            if pdraw < paccept
                #DO accept the proposed values of alpha, beta, gamma
                newalpha = copy(alpha);
                newbeta = copy(beta);
                newgamma = copy(gamma);
            else
                #DO NOT accept the proposed values of alpha, beta, gamma
                #So record the old values prior to proposal
                newalpha = copy(oldalpha);
                newbeta = copy(oldbeta);
                newgamma = copy(oldgamma);
            end
            
            #save parameter set
            global alphavec = push!(alphavec,newalpha);
            global betavec = push!(betavec,newbeta);
            global gammavec = push!(gammavec,newgamma);
            
            #Change the temperature
            # deltaL = (L[end] - L[end-1])^2;
            # if isinf(deltaL) == false
            #     temp = minimum([temp*(deltaL),100.]);
            # end
            temp = 100/(1 + log(1 + 0.1*tic));
            
            global tvec = push!(tvec,temp);
            
            tic += 1;
        end #end while
    end # end let
R"plot($L)"
[alphavec[end],betavec[end],gammavec[end]]
            
            # global cumnum = 0.0;
            # global cumdenom = 0.0;
            # cnumerator = Array{Float64}(undef,1);
            # cdenominator = Array{Float64}(undef,1);
            # let cumnum = 0.0, cumdenom = 0.0
            #     for i=1:numpred
            #         for j=1:numprey
            # 
            #             aij = copy(A[j,i]);
            #             mi = massvec_pred[i]; #i is predator
            #             mj = massvec_prey[j]; #j is prey
            # 
            #             pij = exp(alphamax + betamax*log(mi/mj) + gammamax*log(mi/mj)^2)/(1+exp(alphamax + betamax*log(mi/mj) + gammamax*log(mi/mj)^2));
            # 
            #             num = aij*pij + (1-aij)*(1-pij);
            #             denom = aij + (1-aij);
            # 
            # 
            #             cumnum += num;
            #             cumdenom += denom;
            #         end
            #     end
            #     cnumerator[1] = copy(cumnum);
            #     cdenominator[1] = copy(cumdenom);
            # end
            # 
            # fc_new = cnumerator[1]/cdenominator[1];
            # 
            # #evaluate change in fc
            # deltafc = (fc_old - fc_new)^2;
            
    
    
    
    return L, thetavec, fc
    
end

                