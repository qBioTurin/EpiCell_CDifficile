
cd '/home'/riccardo/cobratoolbox/

initCobraToolbox

cd '/home/riccardo/Documents'/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/

CD196 = importdata('CD196.mat');

% the following steps are aimed at generate reaction/s, and then add the reaction/s to the model. 
% We'll use the 'sink_pheme[c]' reaction from the CD196 model using "sink_gthrd(c)" as demo reaction

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
BIGGdataCD196react = readtable("/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196/BIGGdata_CD196_react.tsv", opts);
clear opts

% demo reactions

data_gthrd = BIGGdataCD196react(strcmp(BIGGdataCD196react.abbreviation, 'sink_gthrd(c)'), :);

cd '/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196'

% create the model with new reactions
% sink_pheme(c)

[CD196SinkHeme] = addReactionGEM_Unito(CD196, {'sink_pheme(c)'}, ...
{'sink reaction for heme'}, ...
cellstr('pheme[c] <==> '), ...
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

cd '/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink'

save('CD196HemeSink.mat', 'CD196SinkHeme');