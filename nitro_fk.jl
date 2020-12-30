using Distributed
using SharedArrays

@everywhere using ProgressMeter
@everywhere using LightGraphs
@everywhere using Distributions
@everywhere using RCall
# @everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nichemodelweb.jl");
# @everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/quantitativeweb.jl");
# @everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nitromax.jl");

@everywhere include("$(homedir())/oceansfoodweb/src/nichemodelweb.jl");
@everywhere include("$(homedir())/oceansfoodweb/src/quantitativeweb.jl");
@everywhere include("$(homedir())/oceansfoodweb/src/nitromax.jl");

@everywhere function calcerror(Q,est_tl)
    # q = UnipartiteQuantiNetwork(Q);
    # pred_tl = trophic_level(q);
    R"pred_tl <- TrophInd($(transpose(Q)))[,1]"
    @rget pred_tl;
    # err = mean(sqrt.((pred_tl .- est_tl).^2));
    err = sum((pred_tl .- est_tl).^2)
    return pred_tl, err
end



reps = 50;
numSp=50;
C=0.1;
minlink = 0.25
steps = 50000;
fkvec = collect(0:0.1:1);

parametervec = [repeat(collect(1:length(fkvec)),inner=reps) repeat(collect(1:reps),outer=length(fkvec))];
its = size(parametervec)[1];

sr2vec = SharedArray{Float64}(length(fkvec),reps);
er2vec = SharedArray{Float64}(length(fkvec),reps);
sr2vecukn = SharedArray{Float64}(length(fkvec),reps);
er2vecukn = SharedArray{Float64}(length(fkvec),reps);

# @showprogress 1 "Computing..." 
@sync @distributed for ii=1:its

    i = parametervec[ii,1];
    r = parametervec[ii,2];
    
    fk = fkvec[i];

    S,
    links,
    unknownlinks,
    Q,
    startQ,
    endQ,
    obs_tl,
    start_tl,
    end_tl,
    errvec,
    Qerrvec,
    tempvec = nitromax(numSp,C,steps,fk,minlink);

    #measure just the stronger links (0.4 to 1)
    #i.e. we are measuring how much we are improving (hopefully) just links with true weights > 0.4
    slinks = findall(x->x>0.4,Q[unknownlinks]);

    R"""
    startr2 <- summary(lm($(startQ[links][slinks]) ~ $(Q[links][slinks])))$adj.r.squared
    endr2 <- summary(lm($(endQ[links][slinks]) ~ $(Q[links][slinks])))$adj.r.squared
    startr2ukn <- summary(lm($(startQ[unknownlinks][slinks]) ~ $(Q[unknownlinks][slinks])))$adj.r.squared
    endr2ukn <- summary(lm($(endQ[unknownlinks][slinks]) ~ $(Q[unknownlinks][slinks])))$adj.r.squared
    """
    @rget startr2;
    @rget endr2;
    @rget startr2ukn;
    @rget endr2ukn;

    sr2vec[i,r] = startr2;
    er2vec[i,r] = endr2;
    sr2vecukn[i,r] = startr2ukn;
    er2vecukn[i,r] = endr2ukn;

end

mer2 = vec(mean(er2vec,dims=2));
msr2 = vec(mean(sr2vec,dims=2));

mer2ukn = vec(mean(er2vecukn,dims=2));
msr2ukn = vec(mean(sr2vecukn,dims=2));


R"""
par(mfrow=c(2,2))
plot($obs_tl,$(start_tl),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($obs_tl,$(end_tl),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($(Q[links]),$(startQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5)
plot($(Q[links]),$(endQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5)
"""

R"""
par(mfrow=c(1,2))
plot($errvec,xlab='Annealing time',ylab='Trophic level error',type='l',lwd=0.5)
plot($Qerrvec,xlab='Annealing time',ylab='Interaction strength error',type='l',lwd=0.5)
"""


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
