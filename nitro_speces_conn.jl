using Distributed
using DataFrames
using Optim
using ShowProgress

@everywhere using LightGraphs
@everywhere using EcologicalNetworks
@everywhere using Distributions
@everywhere using RCall
@everywhere using SharedArrays
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nichemodelweb.jl");
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/quantitativeweb.jl");
@everywhere include("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/src/nitromax.jl");

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

steps = 10000;

spvec = collect(20:1:30);
cvec = collect(0.05:0.01:0.1);
reps = 25;

lspvec = length(spvec);
lcvec = length(cvec);


parametervec1 = [repeat(collect(1:lspvec),outer=lcvec) repeat(collect(1:lcvec),inner=lspvec)];
parametervec = [repeat(parametervec1,reps) repeat(collect(1:reps),inner=lcvec*lspvec)];
its = size(parametervec)[1];

sr2vec = SharedArray{Float64}(lspvec,lcvec,reps);
er2vec = SharedArray{Float64}(lspvec,lcvec,reps);


@showprogress 1 "Computing..." @sync @distributed for ii=1:its
    i = parametervec[ii,1];
    j = parametervec[ii,2];
    r = parametervec[ii,3];
    numSp = Int64(spvec[i]);
    C = cvec[j];
    # r = parametervec[i,3];

    S,
    links,
    Q,
    startQ,
    endQ,
    obs_tl,
    start_tl,
    end_tl,
    errvec,
    Qerrvec,
    tempvec = nitromax(numSp,C,steps);

    R"""
    startr2 <- summary(lm($(startQ[links]) ~ $(Q[links])))$adj.r.squared
    endr2 <- summary(lm($(endQ[links]) ~ $(Q[links])))$adj.r.squared
    """
    @rget startr2;
    @rget endr2;

    sr2vec[i,j,r] = startr2;
    er2vec[i,j,r] = endr2;

end

sr2array = mean(sr2vec,dims=3)[:,:,1];
er2array = mean(er2vec,dims=3)[:,:,1];

R"image($(er2array)-$(sr2array))"

R"hist(($(er2array)/$(sr2array)))"

R"plot($(vec(sr2array)),$(vec(er2array)),pch='.')
