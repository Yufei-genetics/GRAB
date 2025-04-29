#' RPGWAS method in GRAB package
#' 
#' RPGWAS method is an empirical approach to jointly analyzing a binary trait and its risk prediction for related samples in a large-scale biobank.
#' 
#' @details 
#' Additional list of \code{control} in \code{RPGWAS.NullModel()} function.
#' 
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#' 
#' @examples
#' # Step 1: fit a null model
#' PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData = data.table::fread(PhenoFile, header = T)
#' obj.RPGWAS = RPGWAS.NullModel(SparseGRMFile = SparseGRMFile,
#'                               PairwiseIBDFile = PairwiseIBDFile,
#'                               LambdaGenoFile = LambdaGenoFile,
#'                               control = list(ControlOutlier = FALSE))
#'
#' # Step 2: perform score tesinstall.packages("cli")
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/RPGWASMarkers.txt")
#' GRAB.Marker(objNull = obj.RPGWAS,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#' head(read.table(OutputFile, header = T))
#' @export
GRAB.RPGWAS = function(){
  print("Check ?GRAB.RPGWAS for more details about 'RPGWAS' method.")
}

################### This file includes the following functions

# ------------- This file includes the following functions
# 1. checkControl.Marker.RPGWAS(control)
# 2. setMarker.RPGWAS(objNull, control)
# 3. mainMarker.RPGWAS()

# check the control list in marker-level testing
checkControl.Marker.RPGWAS = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         zeta = 0,
                         tol = 1e-5)
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  return(control)
}

checkControl.RPGWAS.NullModel = function(control,
                                         Pheno,
                                         GenoLambda,
                                         SparseGRM,
                                         PairwiseIBD)
{
  default.control = list(MaxQuantile = 0.75,
                         MinQuantile = 0.25,
                         OutlierRatio = 1.5,
                         ControlOutlier = TRUE,
                         MaxNuminFam = 5,
                         MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
                         CovarColList = NULL)
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  if(control$MaxQuantile < control$MinQuantile)
    stop("MaxQuantile(default is 0.75) should be larger than MinQuantile(default is 0.25).")
  
  if(control$OutlierRatio < 0)
    stop("OutlierRatio should be larger than or equal 0 (default is 1.5).")
  
  if(!all(c("SubjID", "Tar", "Risk") %in% colnames(Pheno)))
    stop("The column names of Pheno should at least include ['SubjID', 'Tar', 'Risk'].")
  
  if(any(colnames(SparseGRM) != c("ID1", "ID2", "Value")))
    stop("The column names of SparseGRM should be ['ID1', 'ID2', 'Value'].")
  
  if(any(colnames(PairwiseIBD) != c("ID1", "ID2", "pa", "pb", "pc")))
    stop("The column names of PairwiseIBD should be ['ID1', 'ID2', 'pa', 'pb', 'pc'].")
  
  if(is.null(rownames(GenoLambda)))
    stop("GenoMat for calculating lambda should have rownames.")
  
  if(ncol(GenoLambda) < 100)
    stop("At least 100 SNPs are required for estimating lambda.")
  
  if(length(which(is.na(GenoLambda))) != 0)
    stop("There exists NA values in GenoMat for estimating lambda.")
   
  SubjID.In.Pheno = Pheno$SubjID
  SubjID.In.GenoLambda = rownames(GenoLambda)
  SubjID.In.GRM = unique(c(SparseGRM$ID1, SparseGRM$ID2))
  SubjID.In.IBD = unique(c(PairwiseIBD$ID1, PairwiseIBD$ID2))
  
  if(any(!SubjID.In.Pheno %in% SubjID.In.GenoLambda))
    stop("At least one subject in phenotype matrix does not have GenoLambda information.")
    
  if(any(!SubjID.In.Pheno %in% SubjID.In.GRM))
    stop("At least one subject in phenotype matrix does not have GRM information.")
  
  if(any(!SubjID.In.IBD %in% SubjID.In.GRM))
    stop("At least one subject has IBD information but does not have GRM information.")
  
  return(control)
}

setMarker.RPGWAS = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setRPGWASobjInCPP(objNull$Tarvec,
                   objNull$Riskvec,
                   objNull$designMat,
                   objNull$GRM,
                   objNull$Resid,
                   objNull$lambda,
                   objNull$gammas,
                   objNull$gamma_riskVec,
                   objNull$inv_tX_X,
                   objNull$inv_tX_X_tX,
                   objNull$t0,
                   objNull$Resid.unrelated.outliers,
                   objNull$sum_R_nonOutlier,
                   objNull$R_GRM_R_nonOutlier,
                   objNull$R_GRM_R_TwoSubjOutlier,
                   objNull$R_GRM_R,
                   objNull$MAF_interval,
                   objNull$TwoSubj_list,
                   objNull$ThreeSubj_list,
                   control$SPA_Cutoff,
                   control$zeta,
                   control$tol)
  
  print(paste0("The current control$nMarkersEachChunk is", control$nMarkersEachChunk, ".")) # This file is in 'control.R'
}

mainMarker.RPGWAS = function(genoType, genoIndex, outputColumns)
{
  OutList = mainMarkerInCPP("RPGWAS", genoType, genoIndex);
  
  # print(lapply(OutList, length))
  
  obj.mainMarker = data.frame(Marker = OutList$markerVec,            # marker IDs
                              Info = OutList$infoVec,                # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,          # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,      # alternative allele counts
                              MissingRate = OutList$missingRateVec,  # missing rate
                              zScore = OutList$zScore,               # standardized score statistics
                              Pvalue = OutList$pvalVec,              # marker-level p-value
                              hwepval = OutList$hwepvalVec)
  
  return(obj.mainMarker)
}


RPGWAS.NullModel = function(PhenoFile,          # at least four columns: column 1 is subjID, column 2 is Target phenotype, column 3 is Risk prediction, column 4 covariates
                            CovaColname,
                            PlinkFile,
                            SparseGRMFile,
                            PairwiseIBDFile,
                            control = list(MaxQuantile = 0.75,
                                           MinQuantile = 0.25,
                                           OutlierRatio = 1.5,
                                           MaxNuminFan = 5,
                                           MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)))
{
  Pheno = data.table::fread(PhenoFile)
  
  GenoList = GRAB.ReadGeno(PlinkFile, control = list(AllMarkers = TRUE, ImputeMethod = 'mean'))
  GenoLambda = GenoList$GenoMat
  
  SparseGRM = data.table::fread(SparseGRMFile)
  PairwiseIBD = data.table::fread(PairwiseIBDFile)
  
  Pheno$SubjID = as.character(Pheno$SubjID)
  SparseGRM$ID1 = as.character(SparseGRM$ID1); SparseGRM$ID2 = as.character(SparseGRM$ID2)
  PairwiseIBD$ID1 = as.character(PairwiseIBD$ID1); PairwiseIBD$ID2 = as.character(PairwiseIBD$ID2)
  
  control = checkControl.RPGWAS.NullModel(control, Pheno, GenoLambda, SparseGRM, PairwiseIBD)
  
  MaxQuantile = control$MaxQuantile
  MinQuantile = control$MinQuantile
  OutlierRatio = control$OutlierRatio
  ControlOutlier = control$ControlOutlier
  MaxNuminFam = control$MaxNuminFam
  MAF_interval = control$MAF_interval
  
  SubjID = Pheno$SubjID
  GenoLambda = GenoLambda[match(SubjID, rownames(GenoLambda)),]
  SparseGRM = SparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  PairwiseIBD = PairwiseIBD %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  
  # calculate residual information
  if(is.null(CovaColname)){
    designMat = Pheno %>% select(setdiff(colnames(Pheno), c("SubjID", "Tar", "Risk")))
    designMat = as.matrix(designMat)
  } else{
    designMat = Pheno %>% select(all_of(CovaColname))
    designMat = as.matrix(designMat)
  }
  
  TarVec = Pheno$Tar
  RiskVec = Pheno$Risk
  
  null.fitter = calculate_nullR(tarVec = TarVec,
                                riskVec = RiskVec,
                                covMat = designMat,
                                GenoLambda = GenoLambda)
  
  # Use residual information to define outliers / non-outliers
  ResidMat = null.fitter$ResidMat
  Resid = ResidMat$Resid
  Quant = quantile(Resid, probs = c(MinQuantile, MaxQuantile))
  Range = max(Quant) - min(Quant)
  cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
  
  cat("cutoffVec:\t",cutoffVec,"\n")
  ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2], TRUE, FALSE)
  
  if(ControlOutlier)
  {
    cat("By default, ControlOutlier = TRUE to keep the outliers within a certain range to improve computational efficiency.\n")
    
    while(sum(ResidMat$Outlier) == 0)
    {
      OutlierRatio = OutlierRatio * 0.8
      cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      cat("cutoffVec:\t",cutoffVec,"\n")
      ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                TRUE, FALSE)
      cat("The number of outlier is:", sum(ResidMat$Outlier),"\n")
    }
    
    while(sum(ResidMat$Outlier)/nrow(ResidMat) > 0.05)
    {
      OutlierRatio = OutlierRatio + 0.5
      cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      cat("cutoffVec:\t",cutoffVec,"\n")
      ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                TRUE, FALSE)
      cat("The number of outlier is:", sum(ResidMat$Outlier),"\n")
    }
  }
  
  cat("Outliers information is as below\n")
  print(ResidMat %>% filter(Outlier == TRUE) %>% dplyr::select(SubjID, Resid, Outlier) %>% arrange(Resid))
  
  # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
  SparseGRM1 = SparseGRM
  SparseGRM1$indice1 = match(SparseGRM$ID1, ResidMat$SubjID)
  SparseGRM1$indice2 = match(SparseGRM$ID2, ResidMat$SubjID)
  
  SparseGRM1$pos1 = ResidMat$Resid[SparseGRM1$indice1]
  SparseGRM1$pos2 = ResidMat$Resid[SparseGRM1$indice2]
  SparseGRM1 = SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))
  
  edges = t(SparseGRM1[, c("ID1", "ID2")])
  graph_GRM = make_graph(edges, directed = F)
  graph_list_all = graph_GRM %>% decompose()
  graph_length = lapply(graph_list_all, length)
  
  graph_list_1 = graph_list_all[graph_length == 1]
  SubjID.unrelated = lapply(graph_list_1, get.vertex.attribute) %>% unlist(use.names = FALSE)
  ResidMat.unrelated = ResidMat %>% filter(SubjID %in% SubjID.unrelated)
  SubjID.unrelated.nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(SubjID) %>% unlist(use.names = F)
  
  # Values used in association analyses
  gammas = null.fitter$fit.null[grep("^gamma", names(null.fitter$fit.null))]
  
  gamma.riskVec = gammas[names(gammas) == "gamma.riskVec"]
  gammas = gammas[names(gammas) != "gamma.riskVec"] 
    
  X = cbind(1, designMat)
  tX_X = t(X) %*% X
  inv_tX_X = solve(tX_X)
  inv_tX_X_tX = inv_tX_X %*% t(X)
  
  R_GRM_R = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated) %>% select(Cov) %>% sum
  sum_R_nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(Resid) %>% sum
  R_GRM_R_nonOutlier = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(Cov) %>% sum
  Resid.unrelated.outliers = ResidMat.unrelated %>% filter(Outlier == TRUE) %>% select(Resid) %>% unlist(use.names = F)
  R_GRM_R_TwoSubjOutlier = 0; TwoSubj_list = ThreeSubj_list = list();
  
  # initialize parameters
  graph_list_updated = list()
  graph_list = graph_list_all[graph_length > 1]
  nGraph = length(graph_list)
  index.outlier = 1
  
  if(nGraph != 0)
  {
    cat("Start process the related residual information.\n")
    
    for(i in 1:nGraph)
    {
      if(i %% 1000 == 0)
        cat("Processing the related residual information:\t", i,"/",nGraph,"\n")
      
      comp1 = graph_list[[i]]
      comp3 = V(comp1)$name
      
      # Step 0: calculate variance for the family
      pos1 = match(comp3, SubjID)
      outlierInFam = any(ResidMat$Outlier[pos1])
      
      block_GRM = make.block.GRM(comp1, SparseGRM)
      
      R_GRM_R.temp = as.numeric(t(Resid[pos1]) %*% block_GRM %*% Resid[pos1])
      R_GRM_R = R_GRM_R + R_GRM_R.temp
      
      if(!outlierInFam)
      {
        sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos1])
        R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        next
      }
      
      # cat("Family ", i, " (with outliers) includes ", length(comp3), " subjects:", comp3, "\n")
      
      vcount = vcount(comp1)   # number of vertices 
      
      if(vcount <= MaxNuminFam)
      {
        graph_list_updated[[index.outlier]] = comp1
        index.outlier = index.outlier + 1
        next
      }
      
      # Step 1: remove the edges until the largest family size is <= MaxNuminFam, default is 5.
      
      comp1.temp = comp1
      tempGRM1 = SparseGRM1 %>% filter(ID1 %in% comp3 | ID2 %in% comp3) %>% arrange(Cov)
      for(j in 1:nrow(tempGRM1))
      {
        # cat("j:\t",j,"\n")
        edgesToRemove = paste0(tempGRM1$ID1[j],"|",tempGRM1$ID2[j])
        comp1.temp = delete.edges(comp1.temp, edgesToRemove)
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        if(max(vcount) <= MaxNuminFam)
          break;
      }
      
      # cat("Edge removal complete. Counts of vertices:\t", vcount,"\n")
      
      # Step 2: add the (removed) edges while keeping the largest family size <= MaxNuminFam, default is 5.
      
      tempGRM1 = tempGRM1[1:j,] %>% arrange(desc(Cov))
      comp1 = comp1.temp
      for(k in 1:nrow(tempGRM1))
      {
        # cat("k:\t",k,"\n")
        edgesToAdd = c(tempGRM1$ID1[k], tempGRM1$ID2[k])
        comp1.temp = add.edges(comp1, edgesToAdd)
        
        vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
        # cat("vcount:\t",vcount,"\n")
        
        if(max(vcount) <= MaxNuminFam)
          comp1 = comp1.temp
      }
      
      comp1 = decompose(comp1)
      
      # cat("Edge add complete. Counts of vertices:\t", comp1 %>% sapply(vcount),"\n")
      
      for(k in 1:length(comp1))
      {
        comp11 = comp1[[k]]
        comp13 = V(comp11)$name
        
        pos2 = match(comp13, SubjID)
        outlierInFam = any(ResidMat$Outlier[pos2])
        
        block_GRM = make.block.GRM(comp11, SparseGRM)
        
        R_GRM_R.temp = as.numeric(t(Resid[pos2]) %*% block_GRM %*% Resid[pos2])
        
        if(!outlierInFam){
          sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos2])
          R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
        }else{
          graph_list_updated[[index.outlier]] = comp11
          index.outlier = index.outlier + 1;
        }
      }
    }
    
    cat("Start process the Chow-Liu tree.\n")
    
    # Make a list of array index.
    arr.index = list()
    for(n in 1:MaxNuminFam)
    {
      temp = c()
      for(i in 1:n)
      {
        indexString = rep("c(1, 1, 1)", n)
        indexString[i] = "0:2"
        indexString = paste0(indexString, collapse = "%o%")
        cmd = paste0("temp = c(temp, list(arr.index", i, "=", indexString, "))")
        eval(parse(text = cmd))
      }
      arr.index[[n]] = temp
    }
    
    # build chou-liu-tree.
    n.outliers = length(graph_list_updated)
    if(n.outliers != 0)
    {
      ## The below values are only used in chou.liu.tree
      TwofamID.index = ThreefamID.index = 0
      for(index.outlier in 1:n.outliers)
      {
        if(index.outlier %% 1000 == 0)
          cat("Processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
        
        comp1 = graph_list_updated[[index.outlier]]
        comp3 = V(comp1)$name
        n1 = length(comp3)
        pos3 = match(comp3, SubjID)
        
        Resid.temp = ResidMat$Resid[pos3]
        
        if(n1 == 1)
        {
          Resid.unrelated.outliers = c(Resid.unrelated.outliers, Resid.temp)
          next;
        }
        
        block_GRM = make.block.GRM(comp1, SparseGRM)
        
        tempIBD = PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)
        
        if(n1 == 2)
        {
          TwofamID.index = TwofamID.index + 1;
          
          R_GRM_R_TwoSubjOutlier.temp = as.numeric(t(Resid.temp) %*% block_GRM %*% Resid.temp)
          R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp
          
          Rho.temp = tempIBD$pa + 0.5*tempIBD$pb
          midterm = sqrt(Rho.temp^2 - tempIBD$pa)
          
          TwoSubj_list[[TwofamID.index]] = list(Resid = Resid.temp,
                                                Rho = c(Rho.temp + midterm, Rho.temp - midterm))
          next;
        }
        
        ThreefamID.index = ThreefamID.index + 1;
        
        CLT = chow.liu.tree(N = n1,
                            IBD = tempIBD,
                            IDs = comp3,
                            MAF_interval = MAF_interval)
        
        stand.S.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid.temp)), rep(3, n1))
        
        ThreeSubj_list[[ThreefamID.index]] = list(CLT = CLT,
                                                  stand.S = c(stand.S.temp))
      }
      cat("Completed processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
    }
  }
  
  obj = list(Tarvec = TarVec, Riskvec = RiskVec, designMat = designMat, GRM = SparseGRM1, Resid = Resid, 
             subjData = SubjID, N = length(SubjID), lambda = null.fitter$lambda, gammas = gammas, gamma_riskVec = gamma.riskVec,
             inv_tX_X = inv_tX_X, inv_tX_X_tX = inv_tX_X_tX, t0 = null.fitter$t0, Resid.unrelated.outliers = Resid.unrelated.outliers,
             R_GRM_R = R_GRM_R, R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier,
             sum_R_nonOutlier = sum_R_nonOutlier, R_GRM_R_nonOutlier = R_GRM_R_nonOutlier,
             TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list, MAF_interval = MAF_interval)
  
  class(obj) = 'RPGWAS_NULL_Model'
  return(obj)
}

calculate_nullR = function(tarVec,
                           riskVec,
                           covMat,
                           GenoLambda){
  # Step 1: Fit null model.
  fit.null = fit_null(tarVec = tarVec,
                      riskVec = riskVec,
                      covMat = covMat)
  
  # Step 2: Lambda
  SubjID = rownames(GenoLambda)
  gammas = fit.null[grep("^gamma", names(fit.null))]
  beta1s = fit.null[grep("^beta", names(fit.null))]
  
  Z = cbind(rep(1, nrow(covMat)), covMat)

  info = data.table::data.table(SubjID = SubjID, Resid = 0)
  
  obs.Y = which(!is.na(tarVec))
  obs.info = info[obs.Y, ]
  miss.info = info[-obs.Y, ]
    
  t0 = Z[obs.Y,] %*% gammas[names(gammas) != "gamma.riskVec"]
  t =  t0 + riskVec[obs.Y] * gammas[names(gammas) == "gamma.riskVec"]
  
  obs.info = obs.info %>% mutate(t = c(t),
                                 f = dnorm(t),
                                 F = pnorm(t),
                                 deno = F * (1 - F))
  
  obs.info$deno[which(obs.info$F == 1)] = 1 * pnorm(-t[which(obs.info$F == 1)]) + 1e-200
  obs.info$deno[which(obs.info$F == 0)] = 1e-200
  
  lambdas = c()
  for (i in 1:ncol(GenoLambda)){
    geno = GenoLambda[, i]
    ia2b2 = - fit.null['gamma.riskVec'] * sum(obs.info$f^2 * geno[obs.Y]^2 / obs.info$deno)
    ib2b2 = sum(geno^2)/fit.null['sigma2_sq'] - fit.null['gamma.riskVec'] * ia2b2
    lambda = ia2b2/ib2b2
    lambdas = c(lambdas, lambda)
    
    if (i %% 100 == 0){
      cat('Have finished', i, '/ 1000 snps.\n')
    }
    # print(lambda)
  }
  est_lambda = mean(lambdas)
  
  # Step2: calculate residuals.
  obs.info = obs.info %>% mutate(R_alpha2 = f * (tarVec[obs.Y] - F) / deno,
                                 R2_beta2 = - gammas['gamma.riskVec'] * R_alpha2,
                                 Resid = R_alpha2 - est_lambda * R2_beta2)
  
  
  R1_beta2 = (1/fit.null['sigma2_sq']) * (riskVec - Z %*% beta1s)

  resid.info = rbind(obs.info %>% select(SubjID, Resid), miss.info)
  resid.info$SubjID = factor(resid.info$SubjID, levels = info$SubjID)
  resid.info = resid.info %>% arrange(SubjID)
  
  resid_risk = riskVec - Z %*% beta1s
  resid.info$Resid = resid.info$Resid - est_lambda * (1/fit.null['sigma2_sq']) * resid_risk
  
  return(list(ResidMat = resid.info %>% mutate(SubjID, Resid), 
              lambda = est_lambda, 
              fit.null = fit.null, 
              resid_risk = resid_risk,
              t0 = t0))
}

fit_null = function(tarVec,
                    riskVec,
                    covMat){
  # Stage 1 regression.
  fit_1 = lm(riskVec ~ covMat)
  
  coef_1 = summary(fit_1)$coefficients
  
  if (nrow(coef_1) != ncol(covMat) + 1){
    stop("Please check covMat of there exists multicollinearity or one or more covariates are the identity vector (all elements are 1).")
  }
  
  betas = coef_1[, 'Estimate']
  sigma2 = summary(fit_1)$sigma
  
  # Stage 2 regression.
  fit_2 = glm(tarVec ~ covMat + riskVec, family = binomial(link = "probit"), na.action = na.omit)
  
  coef_2 = summary(fit_2)$coefficients
  
  gammas = coef_2[rownames(coef_2) != 'riskVec', 'Estimate']
  gamma.riskVec = coef_2['riskVec', 'Estimate']
  
  names(betas) = paste("beta", names(betas), sep = ".")
  names(gammas) = paste("gamma", names(gammas), sep = '.')
  names(gamma.riskVec) = 'gamma.riskVec'
  
  out = c(betas, sigma2_sq = sigma2^2, gammas, gamma.riskVec)
  return(out)
}

make.block.GRM = function(graph, 
                          GRM)    # three columns: "ID1", "ID2", and "Value"
{
  comp2 = get.data.frame(graph)
  
  # igraph gives an unexpected additional loop, which may change the block GRM
  # the below is to remove the additional loop
  comp2 = comp2[!duplicated(comp2),]
  
  comp3 = V(graph)$name
  
  colnames(GRM) = c("to", "from", "Value")
  
  n1 = nrow(comp2)
  comp2 = merge(comp2, GRM)
  n2 = nrow(comp2)
  
  if(n1 != n2)
    stop("Ask Wenjian Bi (wenjianb@pku.edu.cn) to check why 'n1 != n2'.")
  
  block_GRM = sparseMatrix(i = match(comp2$from, comp3),
                           j = match(comp2$to, comp3),
                           x = comp2$Value,
                           symmetric = T)
  return(block_GRM)
}

chow.liu.tree = function(N,
                         IBD,
                         IDs,
                         MAF_interval)
{
  CLT = c()
  
  # MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  for(index in 1:length(MAF_interval))
  {
    mu = MAF_interval[index]
    
    # p = c(G0, G1, G2)
    p0 = c((1-mu)^2, 2*mu*(1-mu), mu^2)
    
    # p = c(G00, G10, G20, G01, G11, G21, G02, G12, G22)
    pa.allele2 = c((1-mu)^2, 0, 0, 0, 2*mu*(1-mu), 0, 0, 0, mu^2)
    
    pb.allele1 = c((1-mu)^3, mu*(1-mu)^2, 0, mu*(1-mu)^2, mu*(1-mu), mu^2*(1-mu), 0, mu^2*(1-mu), mu^3)
    
    pc.allele0 = c((1-mu)^4, 2*mu*(1-mu)^3, mu^2*(1-mu)^2, 2*mu*(1-mu)^3, 4*mu^2*(1-mu)^2, 
                   2*mu^3*(1-mu), mu^2*(1-mu)^2, 2*mu^3*(1-mu), mu^4)
    
    # calculate entropy I(Gi, Gj). Noting that entropy of unrelated pairs is zero.
    for(j in 1:nrow(IBD))
    {
      pro = IBD$pa[j] * pa.allele2 + IBD$pb[j] * pb.allele1 + IBD$pc[j] * pc.allele0
      
      entropy = sum(pro * log(pro/pc.allele0), na.rm = T)
      IBD$entropy[j] = entropy
    }
    
    # use the "prim" lgorithm to bulid a maximum spanning tree.
    Max_span_tree = IBD %>% graph_from_data_frame(directed = T) %>% 
      mst(weights = - IBD$entropy, algorithm = "prim") %>% get.edgelist() %>% 
      data.table::as.data.table() %>% rename(ID1 = V1, ID2 = V2)
    
    mst.IBD = merge(Max_span_tree, IBD, all.x = T) %>%
      mutate(idxID1 = match(ID1, IDs), idxID2 = match(ID2, IDs))
    
    arr.prob = array(1, dim = rep(3, N))
    for(i in 1:N)
      dimnames(arr.prob)[[i]] = paste0("ID",i,":",0:2) 
    
    vec = c(mst.IBD$idxID1, mst.IBD$idxID2); vec = vec[duplicated(vec)]
    
    for(k in 1:(N - 1))
    {
      pro = mst.IBD$pa[k] * pa.allele2 + mst.IBD$pb[k] * pb.allele1 + mst.IBD$pc[k] * pc.allele0
      
      matrix.prob = matrix(pro, 3, 3)
      matrix.index1 = mst.IBD$idxID1[k]; matrix.index2 = mst.IBD$idxID2[k]
      for(i in 1:3){
        for(j in 1:3){
          indexString = rep("", N)
          indexString[matrix.index1] = i
          indexString[matrix.index2] = j
          indexString = paste0(indexString, collapse = ",")
          cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] * matrix.prob[",i,",",j,"]")
          # "arr.prob[1,1,] = arr.prob[1,1,] * matrix.prob[1,1]"
          eval(parse(text = cmd))
        }
      }
    }
    
    for(k in 1:(N - 2))
    {
      vector.prob = p0
      vector.index = vec[k]
      for(i in 1:3){
        indexString = rep("", N)
        indexString[vector.index] = i
        indexString = paste0(indexString, collapse = ",")
        cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] / vector.prob[",i,"]")
        # arr.prob[,,1] = arr.prob[,,1] / vector.prob[1]"
        eval(parse(text = cmd))
      }
    }
    
    CLT = cbind(CLT, c(arr.prob))
  }
  
  return(CLT)
}
