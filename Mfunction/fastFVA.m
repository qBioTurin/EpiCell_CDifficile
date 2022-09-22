

cd /home/riccardo/cobratoolbox
initCobraToolbox

model = importdata('/home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Input/CDmodels/CD196HemeSink/CD196HemeSink.mat');
model = buildRxnGeneMat(model);
model.rxnNames(cellfun('isempty', model.rxnNames)) = {'NA'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Random model %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelRand = model;
rxnNameList = modelRand.rxns;

lb = -randi(10, length(modelRand.lb), 1);
ub = randi(10, length(modelRand.ub), 1);
modelRand = changeRxnBounds(modelRand, rxnNameList, lb, 'l');
modelRand = changeRxnBounds(modelRand, rxnNameList, ub, 'u');

[modelRandIrrev] = convertToIrreversible(modelRand);

diseaseState = 'Random';
[minFluxRandom, maxFluxRandom, VminRandom, VmaxRandom] = fluxVariability(modelRandIrrev, 90, 'max');

cd /home/riccardo/Documents/ClostridiumDiff_FBAandPN/EpiCell_CDifficile/Results/

save(['minFlux' diseaseState '.mat'],'minFluxRandom','-v7.3');
save(['maxFlux' diseaseState '.mat'],'maxFluxRandom','-v7.3');
save(['fvamin' diseaseState '.mat'],'VminRandom','-v7.3');
save(['fvamax' diseaseState '.mat'],'VmaxRandom','-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% No Drug model %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modelNoDrug = model;
[modelNoDrugIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(modelNoDrug);

modelTherapy= model;
[modelTherapyIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(modelTherapy);

diseaseState = 'NoDrug';
[minFluxNo,maxFluxNo,optsolNo,retNo,fbasolNo,fvaminNo,fvamaxNo,statussolmin,statussolmax] =...
    fastFVA(modelNoDrugIrrev,90,'max');

save(['minFlux' diseaseState '.mat'],'minFluxNo','-v7.3');
save(['maxFlux' diseaseState '.mat'],'maxFluxNo','-v7.3');
save(['fvamin' diseaseState '.mat'],'fvaminNo','-v7.3');
save(['fvamax' diseaseState '.mat'],'fvamaxNo','-v7.3');

diseaseState = 'Therapy';
[minFluxT,maxFluxT,optsolT,retT,fbasolT,fvaminT,fvamaxT,statussolmin,statussolmax] =...
    fastFVA(modelTherapyIrrev,90,'max');

save(['minFlux' diseaseState '.mat'],'minFluxT','-v7.3');
save(['maxFlux' diseaseState '.mat'],'maxFluxT','-v7.3');
save(['fvamin' diseaseState '.mat'],'fvaminT','-v7.3');
save(['fvamax' diseaseState '.mat'],'fvamaxT','-v7.3');

% mavolcanoplot MATLAB function