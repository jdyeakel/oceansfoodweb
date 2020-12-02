using Distributed
using SharedArrays
@everywhere using Distributions
using CSV
using RCall
using DataFrames

data = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/edgelist.csv",header=true);
massdata_raw = CSV.read("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/masslist.csv",header=true);

preds = data[:Predator];
prey = data[:Prey];
splist = unique([preds;prey]);

upreds = unique(preds);
uprey = unique(prey);

#NOTE: turn into Square matrix
A = Array{Int64}(undef,length(splist),length(splist));
A .= 0;

for i=1:length(splist)
    loc_consumer = findall(x->x==splist[i],preds);
    for j=1:length(loc_consumer)
        rowloc = findall(x->x==prey[loc_consumer[j]],splist)[1];
        A[rowloc,i] = 1;
    end
end


#Make dataframe
df = DataFrame(A);
names!(df,Symbol.(splist))
insert!(df,1,splist,:species)
 
CSV.write("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/adjmatrix.csv",df);

#Create mass list
# massvec = repeat([0],length(splist));
# dfmass = DataFrame(species = splist, mass = repeat([0],length(splist)))
# CSV.write("$(homedir())/Dropbox/Funding/20_NSF_OCE/oceanfoodwebs_2020/foodweb_model/foodwebdata_raw/Antarctica_PotterCove/masslist.csv",dfmass);

names!(dfmass,:species);
insert!(dfmass,2,repeat([0],length(splist)),:mass)