using Distributed
using DataFrames
using Optim


@everywhere using LightGraphs
@everywhere using EcologicalNetworks
@everywhere using Distributions
@everywhere using RCall
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nichemodelweb.jl");
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/quantitativeweb.jl");


# S = 20; C = 0.05 gives good results
#Build food web
S=20;
C=0.05;
@time A,n = nichemodelweb(S,C);
S = size(A)[1];


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
R"plot($tl,$obs_tl)"

# #NOTE: THIS PART IS BROKEN NOW
# g = UnipartiteNetwork(A);
# q = UnipartiteQuantitativeNetwork(Q);
# tl = fractional_trophic_level(g);
# #The observed trophic level of species weighted based on assigned prop contr food to diet in Q
# obs_tl = trophic_level(q);


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

println("Rest easy: no loops here...")

# R"""
# par(mfrow=c(1,2))
# plot($obs_tl,$obs_dN)
# plot($obs_tl,$est_tl)
# """

#measure error
@everywhere function calcerror(Q,est_tl)
    # q = UnipartiteQuantiNetwork(Q);
    # pred_tl = trophic_level(q);
    R"pred_tl <- TrophInd($(transpose(Q)))[,1]"
    @rget pred_tl;
    # err = mean(sqrt.((pred_tl .- est_tl).^2));
    err = sum((pred_tl .- est_tl).^2)
    return pred_tl, err
end


# x0 = Q;
# results_bg = optimize(x->calcerror(x,est_tl),x0,NelderMead());



#locations of links per species
linksperspecies = Array{Array}(undef,S);
for i=1:S
    linksperspecies[i] = findall(!iszero,Q[i,:]);
end

#REWRITE THIS USING OPTIMIZE FUNCTION


#Begin simulated annealing algorithm
reps = 50000;
global temperature = 5;
coolingrate = 0.1;
links = findall(!iszero,Q);
tempvec = Array{Float64}(undef,reps);
tlvec = Array{Float64}(undef,S,reps);
errvec =  Array{Float64}(undef,reps);
Qerrvec = Array{Float64}(undef,reps);
#Randomly generate an initial guess
global estQ = quantitativeweb(A);
# estQ = convert(Array{Float64},A);
# estQ[1:length(estQ)] = estQ[1:length(estQ)] .* (repeat(1./convert.(Float64,length.(linksperspecies)),outer=size(Q)[1]));
# estQ[find(isnan,estQ)] = 0;
startQ = copy(estQ);
#Calculate the initial error
global pred_tl, err = calcerror(estQ,est_tl);

tlvec[:,1] = pred_tl;
errvec[1] = err;
tempvec[1] = temperature;
# Qerrvec[1] = mean(sqrt.((Q[links].-estQ[links]).^2));
Qerrvec[1] = sum((Q[links].-estQ[links]).^2);
@time for i=2:reps
    
    if mod(i,100) == 0
        println("iteration = ",i,"; error = ",round(err,digits=3),"; temp = ",round(temperature,digits=3))
    end
    # estlinkstrength = estQ[links];
    #Modify link strengths based on temperature
    newdist = Normal(0,temperature); 
    # newdist = Normal(0,0.1);
    
    #for each species, modify links - 1, and scale the last
    newestQ = copy(estQ);
    
    #Choose random consumer
    consumers = findall(x->x>1,length.(linksperspecies));
    
    #Choose number of species to adjust based on temperature
    # numsp = minimum([length(consumers),maximum([1,Int64(round(20*temperature,digits=0))])]);
    numsp = 1;
    sptoadjust = rand(consumers,numsp);
    
    for j=sptoadjust
        slinks = linksperspecies[j];
        #choose link to alter
        linktochange = rand(slinks);
        newestQ[j,linktochange] = maximum([minimum([0.99,estQ[j,linktochange]*(1+rand(newdist))]),0.001]);
        #rescale
        newestQ[j,:] = newestQ[j,:]./(sum(newestQ[j,:]));
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
        global temperature = maximum([0.00001,temperature * ((err_new)/(err))]);
    
        #Accept new changes
        global estQ = copy(newestQ);
        global err = copy(err_new);
    else
        #roll dice
        # rdraw = rand();
        # if rdraw < exp(-(err_new - err) / temperature )
        #     #Adjust temperature
        #     # global temperature = maximum([0.00001,temperature * ((err_new)/(err))]);
        
        #     #Accept new changes
        #     global estQ = copy(newestQ);
        #     global err = copy(err_new);
        # end
    end
    
    # #Parallelized annealing sampler
    # sprob = exp( (err_old - err_new) / temperature );
    # sdraw = rand();
    # if sdraw < sprob
    #     nproc = 80;
    #     parerr = SharedArray{Float64}(nproc);
    #     parnewestQvec = SharedArray{Float64}(nproc,size(Q)[1]);
    #     parsp = SharedArray{Int64}(nproc);
    #     @sync @parallel for p = 1:nproc
    #         #Sample across processors
    #         newdist = Normal(0,1);
    # 
    #         #for each species, modify links - 1, and scale the last
    #         parnewestQ = copy(estQ);
    # 
    #         #Choose random consumer
    #         consumers = find(x->x>1,length.(linksperspecies));
    #         j = rand(consumers);
    #         slinks = linksperspecies[j];
    #         #choose link to alter
    #         linktochange = rand(slinks);
    #         parnewestQ[j,linktochange] = maximum([minimum([0.99,estQ[j,linktochange]*(1+rand(newdist))]),0.001]);
    #         #rescale
    #         parnewestQ[j,:] = parnewestQ[j,:]./(sum(parnewestQ[j,:]));
    # 
    #         parpred_tl_new,parerr_new = calcerror(parnewestQ,est_tl);
    #         parerr[p] = parerr_new;
    #         parnewestQvec[p,:] = parnewestQ[j,:];
    #         parsp[p] = j;
    #     end
    #     minparerr = findmin(parerr)[2];
    #     if parerr[minparerr] < err
    #         #Accept new changes
    #         err = copy(err_new);
    #         estQ[parsp[minparerr],:] = parnewestQvec[minparerr,:];
    #     end
    # end
    # 
    
    
    #record temperature
    tempvec[i] = temperature;
    tlvec[:,i] = pred_tl_new;
    errvec[i] = err;
    #Does this increase the accuracy with which we predict Q?
    # Qerrvec[i] = mean(sqrt.((Q[links].-estQ[links]).^2));
    Qerrvec[i] = sum((Q[links].-estQ[links]).^2);
end

endQ = copy(estQ);

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/errtemp.pdf";
R"""
pdf($namespace,height = 6, width = 10)
par(mfrow=c(1,2))
plot($errvec)        
plot($tempvec)
dev.off()
"""

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/tlerrQerr.pdf";
R"""
pdf($namespace,height = 6, width = 10)
par(mfrow=c(1,2))
plot($errvec,xlab='Annealing time',ylab='Trophic level error',type='l',lwd=0.5)
plot($Qerrvec,xlab='Annealing time',ylab='Interaction strength error',type='l',lwd=0.5)
dev.off()
"""

namespace = "$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/figures/bacomp.pdf";
R"""
pdf($namespace,height = 10, width = 10)
par(mfrow=c(2,2))
plot($obs_tl,$(tlvec[:,1]),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($obs_tl,$(tlvec[:,reps]),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($(Q[links]),$(startQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5)
plot($(Q[links]),$(endQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5)
dev.off()
"""

R"""
plot($(Q[links]),jitter($(sqrt.((endQ[links] .- Q[links]).^2))))
"""
