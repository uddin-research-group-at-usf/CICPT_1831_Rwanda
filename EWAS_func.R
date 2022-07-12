
##############   get_all_lm(dmrs,indep,df)   ##################
### Gets Accounting for the variables specified in the *indep* parameter, 
### this function runs linear regression for each highly enriched cpgs with the mothers methylation as dependent and her childs methylation as independent


dmr_list <- list()

cpg_lm <- function(ass_cpgs # chr vector of 3
                   ,indep # covariates using lm format in quotes
                   ,df){ 
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    y_var <- as.matrix(y_var)
    
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    x_var <- as.matrix(x_var)
    
    x_var <- as.matrix(x_var)          
    
    lm <- lm(y_var ~ x_var,data = df)
    
    res <- as.data.frame(summary(lm)$coefficients)
    #} 
    
    dmr_list[[i]] <- res
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

fin_list <- list()

get_all_lm <- function(dmrs 
                       ,indep # covariates using lmer format in quotes
                       ,df){
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]    
    
    res <- cpg_lm(ass_cpgs,indep,df)
    
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res
    
  }
  fin_list
}

##############   get_all_rr(dmrs,indep,df,plotdir)   ##################
### This is similar to `get_all_lm()` but implements idge regression for top three enriched CpGs of DMRs in mothers. 

Corner_text <- function(text, location="topright"){
  legend(location,legend=text, bty ="n", pch=NA) 
}

cpg_rr <- function(ass_cpgs # chr vector of 3
                   ,indep # chr string of variables
                   ,df
                   ,plotdir){ 
  
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    y_var <- as.matrix(y_var)
    
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    x_var <- as.matrix(x_var)
    
    lambda_seq <- 10^seq(3, -2, by = -.1)
    
    cv_glmnet <- glmnet::cv.glmnet(x_var,y_var, alpha = 0, lambda  = lambda_seq, nfolds = 5,standardize = FALSE)
    
    best_lambda <- cv_glmnet$lambda.min
    
    best_ridge <- glmnet::glmnet(x_var,y_var, alpha = 0, lambda = best_lambda,standardize = FALSE)
    
    res <- as.data.frame(coef(best_ridge)[,1])
    colnames(res) <- "coef"
    
    
    dmr_list[[i]] <- res
    
    L <- length(cv_glmnet$glmnet.fit$lambda)
    x <- log(cv_glmnet$glmnet.fit$lambda[L])
    y <- cv_glmnet$glmnet.fit$beta[, L]
    labs <- names(cv_glmnet$glmnet.fit$beta[, L])
    
    png(paste0(plotdir,ass_cpgs[i],'.png'))
    plot(cv_glmnet$glmnet.fit, xvar = "lambda", label =TRUE)
    abline(v = log(best_lambda))
    legend("topright", paste(seq(1,9,1),"=",labs ), cex= 0.8)
    title(ass_cpgs[i], line = -2,cex.lab=1.2)
    Corner_text(text=paste("lambda =", best_lambda),location= "bottomright")
    dev.off()
    
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

fin_list <- list()

get_all_rr <- function(dmrs
                       ,indep # covariates using lmer format in quotes
                       ,df
                       ,plotdir){ # specify output directory for trace plots
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]
    
    
    res <- cpg_rr(ass_cpgs,indep,df,plotdir)
    
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res
    
  }
  fin_list
  
}

##############   rr_AIC(dmrs,indep,df)   ##################
### Gets AIC of Ridge Regression model. Arguments are similar to previous functions.

AIC_ridge <- function(ridge_fit){ # Function to get AIC of ridge regression
  tLL <- ridge_fit$nulldev - deviance(ridge_fit)
  k <- ridge_fit$df
  n <- ridge_fit$nobs
  AIC <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  print(AIC)
}


cpg_rr_AIC <- function(ass_cpgs # chr vector of 3
                       ,indep # chr string of variables
                       ,df){ 
  
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    
    
    y_var <- as.matrix(y_var)
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    
    
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    
    x_var <- as.matrix(x_var)
    lambda_seq <- 10^seq(3, -2, by = -.1)
    cv_glmnet <- glmnet::cv.glmnet(x_var,y_var, alpha = 1, lambda  = lambda_seq, nfolds = 5,standardize = FALSE)
    
    best_lambda <- cv_glmnet$lambda.min
    
    best_ridge <- glmnet::glmnet(x_var,y_var, alpha = 1
                                 , lambda = best_lambda,standardize = FALSE)
    
    dmr_list[[i]] <- AIC_ridge(best_ridge)
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

rr_AIC <- function(dmrs, indep, df){
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]
    
    res <- cpg_rr_AIC(ass_cpgs, indep,df)
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res  
  }
  
  fin_list
}


##############   lm_AIC(dmrs,indep,df)   ##################
### Gets AIC of Linear Regression model. Arguments are similar to previous functions.


dmr_list <- list()

cpg_lm_AIC <- function(ass_cpgs 
                       ,indep 
                       ,df){ 
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    y_var <- as.matrix(y_var)
    
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    x_var <- as.matrix(x_var)
    
    x_var <- as.matrix(x_var)          
    
    lm <- lm(y_var ~ x_var,data = df)
    print(AIC(lm))
    
    res <- as.data.frame(AIC(lm))
    
    
    dmr_list[[i]] <- res
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

fin_list <- list()

lm_AIC <- function(dmrs
                   ,indep 
                   ,df){
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]    
    
    res <- cpg_lm_AIC(ass_cpgs,indep,df)
    
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res
    
  }
  fin_list
}

##############   cpg_ss(dmrs,indep,df)   ##################
### Gets Stability Selection Results. Arguments are similar to previous functions.

stabsel_JD <- function(x,error=0.05,type=c("pfer","pcer"),pi_thr=0.6){
  if(pi_thr <= 0.5 | pi_thr >= 1) stop("pi_thr needs to be > 0.5 and < 1!")
  if(class(x$fit)[1]=="multnet"){
    p <- dim(x$fit$beta[[1]])[1]
  }else{
    p <- dim(x$fit$beta)[1]
  }
  type <- match.arg(type)
  switch(type,
         "pcer"={
           if(error>=1 | error<=0)stop("pcer needs to be > 0 and < 1!")
           qv <- ceiling(sqrt(error* p * (2*pi_thr-1)*p)) },
         "pfer"={
           qv <- ceiling(sqrt(error * (2*pi_thr-1)*p)) }
  )
  if(x$qs[length(x$qs)]<=qv){ lpos <- length(x$qs)
  }else{
    lpos <- which(x$qs>qv)[1]
  }
  if(!is.na(lpos)){stable <- which(x$x[,lpos]>=pi_thr)}else{
    stable <- NA
  }
  out <- list(stable=stable,lambda=x$fit$lambda[lpos],lpos=lpos,error=error,type=type)
  return(x$x[,lpos])
}

cpg_ss <- function(ass_cpgs # chr vector of 3
                   ,indep # chr string of variables
                   ,df){ 
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    y_var <- as.matrix(y_var)
    
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    x_var <- as.matrix(x_var)
    
    
    spath <- c060::stabpath(y = y_var, x = x_var , steps = 100,mc.cores = 2,
                            weakness = 1, alpha = 0.001)  
    
    res <- data.frame(stabsel_JD(spath,error = 1 ,type = "pfer", pi_thr = 0.6))
    
    dmr_list[[i]] <- res
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

fin_list <- list()

get_all_ss <- function(dmrs
                       ,indep # covariates using lmer format in quotes
                       ,df){  # specify output directory for trace plots
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]
    
    
    res <- cpg_ss(ass_cpgs,indep,df)
    
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res
    
  }
  fin_list
  
}

##############   rr_emm(dmrs,indep,df)   ##################
### Gets get Ridge EMMs. Arguments are similar to previous functions.

cpg_rr_emm <- function(ass_cpgs # chr vector of 3
                       ,indep # chr string of variables
                       ,df){ 
  
  for(i in seq_along(ass_cpgs)){ 
    
    y_cpg <- paste0("C_",ass_cpgs[i])
    y_var <- scale(df[,names(df) %in% y_cpg])
    y_var <- as.matrix(y_var)
    
    new_dep <- c(paste0("M_",ass_cpgs[i]), indep)
    
    x_var <- as.matrix(df[,names(df) %in% new_dep])
    x_var[,paste0("M_",ass_cpgs[i])] <- scale(x_var[,paste0("M_",ass_cpgs[i])])
    x_var <- as.matrix(x_var)
    
    lambda_seq <- 10^seq(3, -2, by = -.1)
    
    cv_glmnet <- glmnet::cv.glmnet(x_var,y_var, alpha = 0, lambda  = lambda_seq, nfolds = 5,standardize = FALSE)
    
    best_lambda <- cv_glmnet$lambda.min
    
    best_ridge <- glmnet::glmnet(x_var,y_var, alpha = 0, lambda = best_lambda,standardize = FALSE)
    
    
    mm <- data.matrix(x_var, new_df)
    
    pred <- predict(best_ridge,  mm)
    
    mm_df <- as.data.frame(mm)
    
    
    control_rows <- mm_df %>% 
      filter(M_Exposure == 0 )%>%
      rownames()
    
    expos_rows   <- mm_df %>% 
      filter(M_Exposure == 1 )%>%
      rownames()
    
    emm_control <- sum(pred[rownames(pred) %in% control_rows,])/length(control_rows)
    
    emm_exp <- sum(pred[rownames(pred) %in% expos_rows,])/length(expos_rows)
    
    final_res <- data.frame( emm = rbind(emm_control,emm_exp),
                             SE = sigma(best_ridge)/ sqrt(c(length(control_rows),length(expos_rows))) )
    
    
    
    
    
    final_res$group <- c("control","exposed")
    
    dmr_list[[i]] <-  final_res
    
  }
  names(dmr_list) <- ass_cpgs
  dmr_list
}

fin_list <- list()

rr_emm <- function(dmrs
                   ,indep # covariates using lmer format in quotes
                   ,df){ 
  
  new_dmr_list <- setNames(split(dmrs, seq(nrow(dmrs))), rownames(dmrs))
  
  for(i in seq_along(new_dmr_list)){
    
    ass_cpgs <- unlist(strsplit(new_dmr_list[[i]]$leadingEdge," "))   
    ass_cpgs <- substr(ass_cpgs, 1, nchar(ass_cpgs)-1)[1:3]
    
    
    res <- cpg_rr_emm(ass_cpgs,indep,df)
    
    fin_list[[paste(new_dmr_list[[i]]$DMR)]] <- res
    
  }
  fin_list
  
}

##############   combine_emm(rr_emm_output, dmr)   ##################
### Gets get averaged Ridge EMMs between DMRs. Arguments are similar to previous functions.

combine_emm <- function(rr_emm_output, dmr){ # Function get averaged emms between DMRs
  for(i in seq_along(dmr)){
    df <- do.call(rbind,rr_emm_output[[dmr[i]]])
    df$DMR <- names(rr_emm_output[dmr[i]])
    
    exp_df <- data.frame( emm = mean(df[df$group == "exposed",]$emm),
                          se = mean(df[df$group == "exposed",]$SE),
                          group = "exposed",
                          DMR = names(rr_emm_output[dmrs_exps[i]])) 
    cont_df <- data.frame( emm = mean(df[df$group == "control",]$emm),
                           se = mean(df[df$group == "control",]$SE),
                           group = "control",
                           DMR = names(rr_emm_output[dmrs_exps[i]]))
    
    final_df <- rbind(cont_df,exp_df)
    fin_list[[i]] <- final_df 
  }
  fin_com_df <- do.call(rbind, fin_list)
  fin_com_df
} 


##############   save_df(res_list,a = names(res_list))   ##################
### converting list objects from previous functions to df. Assign arguement A for DMR names

save_df <- function(res_list,a = names(res_list)){  
  for(i in seq_along(a)){
    rr_res_save <- do.call(rbind,res_list[[a[[i]]]])
    rr_res_save$var <- sub(".*\\.", "", rownames(rr_res_save))
    rr_res_save$CpG <- sub("\\..*", "", rownames(rr_res_save))
    rr_res_save$DMR <- names(res_list[a[[i]]])
    rownames(rr_res_save) <- NULL
    rr_res_save
    fin_list[[i]] <- rr_res_save
  }
  fin_list_df <- do.call(rbind,fin_list)
  return(fin_list_df)
}
