function nitromax(numSp,C,steps,fk)

    # S = 0;
    # while S <= numSp/2
        A,n = nichemodelweb(numSp,C);
        S = size(A)[1];
    # end
        # println(S)

    
    

    #Set 'true' interaction strengths
    Q = quantitativeweb(A);

    R"""
        library(MASS)  
        library(NetIndices)
        rtl<-TrophInd($(transpose(A)))
        qtl <- TrophInd($(transpose(Q)))
    """
    @rget rtl;
    @rget qtl;
    tl = rtl[:,1];
    obs_tl = qtl[:,1];
    # R"plot($tl,$obs_tl)"


    #Build assigned nitrogen values
    unifdist = function(mu,sigma,S) 
        x = rand(Normal(mu,sigma),S); 
        return(x); 
    end
    #Nitrogen isotope variability
    sigmaN = 0.00001;
    obs_dN = ((obs_tl .- 1) .* 3.5).+(unifdist(0,sigmaN,S));

    #estimated trophic level from the Nitrogen isotopes
    est_tl = (obs_dN ./ 3.5) .+ 1;

    # println("Rest easy: no loops here...")

    # R"""
    # par(mfrow=c(1,2))
    # plot($obs_tl,$obs_dN)
    # plot($obs_tl,$est_tl)
    # """

    # #measure error
    # function calcerror(Q,est_tl)
    #     # q = UnipartiteQuantiNetwork(Q);
    #     # pred_tl = trophic_level(q);
    #     R"pred_tl <- TrophInd($(transpose(Q)))[,1]"
    #     @rget pred_tl;
    #     # err = mean(sqrt.((pred_tl .- est_tl).^2));
    #     err = sum((pred_tl .- est_tl).^2)
    #     return pred_tl, err
    # end


    # x0 = Q;
    # results_bg = optimize(x->calcerror(x,est_tl),x0,NelderMead());



    #locations of links per species
    linksperspecies = Array{Array}(undef,S);
    for i=1:S
        linksperspecies[i] = findall(!iszero,Q[i,:]);
    end

    #Choose random consumer (non-specialist!)
    consumers = findall(x->x>1,length.(linksperspecies));
    # specialistconsumers = findall(x->x==1,length.(linksperspecies));
    
    numberknown = Int64(floor(fk*length(consumers)));
    known = sort(sample(consumers,numberknown,replace=false));

    unknownconsumers = setdiff(consumers,known);
    Qunknown = copy(Q);
    #set known links to zero so that Qunknown is just unknown weights
    Qunknown[known,:] .= 0.;

    #REWRITE THIS USING OPTIMIZE FUNCTION


    #Begin simulated annealing algorithm
    # steps = 50000;
    temperature = 5;
    coolingrate = 0.1;
    links = findall(!iszero,Q);
    unknownlinks = findall(!iszero,Qunknown);
    tempvec = Array{Float64}(undef,steps);
    tlvec = Array{Float64}(undef,S,steps);
    errvec =  Array{Float64}(undef,steps);
    Qerrvec = Array{Float64}(undef,steps);
    #Randomly generate an initial guess
    estQ = quantitativeweb(A);
    # estQ = convert(Array{Float64},A);
    # estQ[1:length(estQ)] = estQ[1:length(estQ)] .* (repeat(1./convert.(Float64,length.(linksperspecies)),outer=size(Q)[1]));
    # estQ[find(isnan,estQ)] = 0;
    
    #Plug in known consumer weights
    if length(known) > 0
        estQ[known,:] .= Q[known,:];
    end

    startQ = copy(estQ);

    #Calculate the initial error
    pred_tl, err = calcerror(estQ,est_tl);

    tlvec[:,1] = pred_tl;
    errvec[1] = err;
    tempvec[1] = temperature;
    # Qerrvec[1] = mean(sqrt.((Q[links].-estQ[links]).^2));
    Qerrvec[1] = sum((Q[links].-estQ[links]).^2);



    for i=2:steps
        
        # if mod(i,1000) == 0
        #     println("iteration = ",i,"; error = ",round(err,digits=3),"; temp = ",round(temperature,digits=3))
        # end

        
        # estlinkstrength = estQ[links];
        #Modify link strengths based on temperature
        newdist = Normal(0,temperature); 
        # newdist = Normal(0,0.1);
        
        #for each species, modify links - 1, and scale the last
        newestQ = copy(estQ);
        
        #Choose random consumer
        # consumers = findall(x->x>1,length.(linksperspecies));
        

        
        # unknowngenconsumers = setdiff(unknownconsumers,specialistconsumers);
        
        #Choose number of species to adjust based on temperature
        # numsp = minimum([length(consumers),maximum([1,Int64(round(20*temperature,digits=0))])]);
        numsp = 1;
        if length(unknownconsumers) > 0
            sptoadjust = sample(unknownconsumers,numsp,replace=false);
            for j=sptoadjust
                slinks = linksperspecies[j];
                #choose link to alter
                linktochange = rand(slinks);
                newestQ[j,linktochange] = maximum([minimum([0.99,estQ[j,linktochange]*(1+rand(newdist))]),0.001]);
                #rescale
                newestQ[j,:] = newestQ[j,:]./(sum(newestQ[j,:]));
            end
        end

        pred_tl_new,err_new = calcerror(newestQ,est_tl);
        
        #Lower the temperature relative to the ratio of err/err_new... if err_new is smaller, temperature is lowered; if err_new is larger, temperature is raised
        
        #SIMULATED ANNEALING
        #Acceptance probability
        # prob = exp( (err - err_new) / temperature );
        # rdraw = rand();
        # if rdraw < prob
        #     #Accept new changes
        #     estQ = copy(newestQ);
        #     err = copy(err_new);
        # end
        # #Adjust temperature
        # # temperature = temperature * (err_new/err);
        # temperature *= 1-coolingrate;
        # 
        
        #METROPOLIS
        #only accept changes if error is lowered
        err_old = copy(err);
        
        if err_new <= err
            #Adjust temperature
             temperature = maximum([0.00001,temperature * ((err_new)/(err))]);
        
            #Accept new changes
            estQ = copy(newestQ);
            err = copy(err_new);
        else
            #roll dice
            # rdraw = rand();
            # if rdraw < exp(-(err_new - err) / temperature )
            #     #Adjust temperature
            #     #  temperature = maximum([0.00001,temperature * ((err_new)/(err))]);
            
            #     #Accept new changes
            #      estQ = copy(newestQ);
            #      err = copy(err_new);
            # end
        end
        
        
        #record temperature
        tempvec[i] = temperature;
        tlvec[:,i] = pred_tl_new;
        errvec[i] = err;
        #Does this increase the accuracy with which we predict Q?
        # Qerrvec[i] = mean(sqrt.((Q[links].-estQ[links]).^2));
        Qerrvec[i] = sum((Q[links].-estQ[links]).^2);
    end

    endQ = copy(estQ);

    start_tl = tlvec[:,1];
    end_tl = tlvec[:,steps];

    return(S,links,unknownlinks,Q,startQ,endQ,obs_tl,start_tl,end_tl,errvec,Qerrvec,tempvec)
end
