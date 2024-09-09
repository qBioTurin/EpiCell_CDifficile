
cd '/home'/riccardo/Documents/Git/cobratoolbox/

initCobraToolbox

cd '/home/riccardo/Documents/FBAandPN/data/CD196'

CD196 = importdata('CD196.mat');

% the following steps are aimed at generate reaction/s, 
% and then add the reaction/s to the model.     

% 1) We’ll use the ‘HEMEti’ reaction from the Cbifer model
% 2) We’ll use the ‘EX_pheme(e)’ reaction from the Cbifer model

cd '/home/riccardo/Documents/FBAandPN/data/Cbifer'

Cbifer = importdata('Cbifer.mat');

opts = delimitedTextImportOptions("NumVariables", 7);
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
% Specify column names and types
opts.VariableNames = ["VarName1", "rxnNames", "abbreviation", "equation", "lowbnd", "uppbnd", "obj_coef"];
opts.VariableTypes = ["double", "char", "char", "char", "double", "double", "double"];
% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
% Specify variable properties
opts = setvaropts(opts, ["rxnNames", "abbreviation", "equation"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["rxnNames", "abbreviation", "equation"], "EmptyFieldRule", "auto");
% Import the data
BIGGdataCbiferreact = readtable("/home/riccardo/Documents/FBAandPN/data/Cbifer/BIGGdata_Cbifer_react.tsv", opts);
clear opts

% demo reactions

dataHEMEti = BIGGdataCbiferreact(strcmp(BIGGdataCbiferreact.abbreviation, 'HEMEti'), :);
dataEX_pheme = BIGGdataCbiferreact(strcmp(BIGGdataCbiferreact.abbreviation, 'EX_pheme(e)'), :);

dataSHEMEabc = BIGGdataCbiferreact(strcmp(BIGGdataCbiferreact.abbreviation, 'SHEMEabc'), :);
dataEX_sheme = BIGGdataCbiferreact(strcmp(BIGGdataCbiferreact.abbreviation, 'EX_sheme(e)'), :);

cd '/home/riccardo/Documents/FBAandPN/data/CD196'

% create the model with new reactions

% HEMEti

[CD196_heme] = addReactionGEM_Unito(CD196, {'HEMEti'}, ...
{'Heme transport via ABC system'}, ...
cellstr('atp[c] + h2o[c] + pheme[e] <==> adp[c] + h[c] + pheme[c] + pi[c]'), ...
1, -1000, 1000, ...
{'Transport, extracellular'}, ...
{'Unknown'}, {'x(21)'}, {''});

% EX_pheme

[CD196_heme] = addReactionGEM_Unito(CD196_heme, {'EX_pheme(e)'}, ...
{'exchange reaction for heme'}, ...
cellstr('pheme[e] <==> '), ...
1, -1000, 1000, ...
{'Exchange/demand reaction'}, ...
{''}, {''}, {''});

% SHEMEabc

[CD196_heme] = addReactionGEM_Unito(CD196_heme, {'SHEMEabc'}, ...
{'Siroheme transport via ABC system'}, ...
cellstr('atp[c] + h2o[c] + sheme[e] <==> adp[c] + h[c] + pi[c] + sheme[c]'), ...
1, -1000, 1000, ...
{'Transport, extracellular'}, ...
{''}, {''}, {''});

% EX_sheme

[CD196_heme] = addReactionGEM_Unito(CD196_heme, {'EX_sheme(e)'}, ...
{'Siroheme exchange'}, ...
cellstr('sheme[e] <==> '), ...
1, -1000, 1000, ...
{'Exchange/demand reaction'}, ...
{''}, {''}, {''});

% [CD196_heme, HTABLE] = ...
%     addReactionGEM(model = CD196, ...
%     rxns = string(dataEX_sheme.abbreviation), ...
%     rxnNames = string(dataEX_sheme.rxnNames), ...
%     rxnFormulas = dataEX_sheme.equation, ...
%     rev = 0, ...
%     lb = dataEX_sheme.lowbnd, ...
%     ub = dataEX_sheme.uppbnd, ...
%    nRxn = dataEX_sheme.VarName1, ...
%    subSystems = Cbifer.subSystems(dataEX_sheme.VarName1), ...
%    grRules = Cbifer.grRules(dataEX_sheme.VarName1), ...
%    rules = Cbifer.rules(dataEX_sheme.VarName1), ...
%    genes = Cbifer.genes(dataEX_sheme.VarName1), HTABLE);

% saving the new model

cd '/home/riccardo/Documents/FBAandPN/data/CD196_heme'

save('CD196_heme.mat', 'CD196_heme');

% loading a model

cd '/home/riccardo/Documents/FBAandPN/data/iCN900'
iCN900 = importdata('iCN900.mat');

spy(iCN900.S);

Sbin = zeros(size(iCN900.S));
Sbin(find(iCN900.S)) = 1;

for i = 1 : length(iCN900.mets)
metConnectivity(i,1) = sum(Sbin(i,:));
end
loglog(sort(metConnectivity,'descend'),'*')
xlabel('metabolite number (rank ordered) - log scale');
ylabel('number of reactions - log scale')
