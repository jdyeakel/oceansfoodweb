using Distributed

@everywhere using ProgressMeter
@everywhere using LightGraphs
@everywhere using Distributions
@everywhere using RCall
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nichemodelweb.jl");
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/quantitativeweb.jl");
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nitromax.jl");
@everywhere function calcerror(Q,est_tl)
    # q = UnipartiteQuantiNetwork(Q);
    # pred_tl = trophic_level(q);
    R"pred_tl <- TrophInd($(transpose(Q)))[,1]"
    @rget pred_tl;
    # err = mean(sqrt.((pred_tl .- est_tl).^2));
    err = sum((pred_tl .- est_tl).^2)
    return pred_tl, err
end


numSp=50;
C=0.05;
steps = 50000;
fk = 0.8;

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
tempvec = nitromax(numSp,C,steps,fk);

R"""
par(mfrow=c(2,2))
plot($obs_tl,$(start_tl),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($obs_tl,$(end_tl),xlab='True trophic level',ylab='Predicted trophic level',pch=16,cex=0.5)
plot($(Q[links]),$(startQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5,col='gray')
points($(Q[unknownlinks]),$(startQ[unknownlinks]),pch=16,cex=0.5,col='black')
plot($(Q[links]),$(endQ[links]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5,col='gray')
points($(Q[unknownlinks]),$(endQ[unknownlinks]),xlab='True interaction strengths',ylab='Predicted interaction strengths',pch=16,cex=0.5,col='black')
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
