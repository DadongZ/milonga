#' Multiple imputation for repeated binary outcomes
#'
#' @param data Data frame or matrix of multivariate binary outcome with missing values to impute, NA denotes missing values
#' @param cols NULL or a character or integer vector giving the columns for imputation
#' @param nprior integer
#' @param fprior.list a list of prior frequencies. See "Details"
#' @param burnin number of burnin iterations, NULL or positive integer
#' @param rep number of iterations for imputation, default is 10
#' @return list of imputed data frames
#' @export
#' @examples
#' data(polio)
#' milonga(polio, cols=2:7)

milonga<-function(data, cols=NULL, nprior=100, fprior.list=NULL, burnin=NULL, rep=10){
  if(is.null(cols)) {
    data<-as.matrix(data)
    } else  {
    data<-as.matrix(data[, cols])
    }
  n<-nrow(data)
  m<-ncol(data)
  pattern<-do.call(expand.grid, replicate(m, 0:1, simplify=FALSE))
  epriors<-getPatternPrior(pattern, fprior.list)
  npriors<-nprior*epriors
  colnames(pattern)<-colnames(data)
  K<-nrow(pattern)
  cdata<-data[complete.cases(data), ]
  
  impi<-function(datai){
      alpha<-NULL
      for(k in 1:K){
            a<-pattern[k,]
            srow<-grepl(paste(a, collapse="-"), 
                        apply(datai, 1, paste, collapse="-"))
                        nkprior<-npriors[k]
            alpha[k]<-sum(srow)+nkprior   
      } 
      ptheta<-rdirichlet(1,alpha)
      datax<-matrix(NA, ncol=m, nrow=n)
      for (i in 1:n){
            nna<-sum(is.na(data[i, ]))
            if(nna==0){
              datax[i, ]<-data[i, ]
            } else if (nna==m){
                  sp<-rmultinom(1,1,ptheta)==1
                  datax[i,]<-as.matrix(pattern[sp,])
            } else {
                  b<-data[i,]
                  b[is.na(b)]<-"."
                  prow<-grepl(paste(b, collapse="-"), 
                              apply(pattern, 1, paste, collapse="-"))
                  sptheta<-ptheta[prow]
                  spattern<-pattern[prow,]
                  ss<-rmultinom(1,1,sptheta/sum(sptheta))==1
                  datax[i,]<-as.matrix(spattern[ss,])
            }
       }
    return(datax)
  }

  if(is.null(burnin) || burnin < 1){
      burninlist=list()
      datalist<-vector("list", rep)
      datalist[[1]]<-impi(cdata)

      for(i in 1:(rep-1)){
          colnames(datalist[[i]])<-colnames(data)
          datalist[[i+1]]<-impi(datalist[[i]])
      }
  }else{
      burninlist<-vector("list", burnin)
      burninlist[[1]]<-impi(cdata)

      for(i in 1:(burnin-1)){
          colnames(burninlist[[i]])<-colnames(data)
          burninlist[[i+1]]<-impi(burninlist[[i]])
      }

      datalist<-vector("list", rep)
      datalist[[1]]<-impi(burninlist[[burnin]])

      for(i in 1:(rep-1)){
          colnames(datalist[[i]])<-colnames(data)
          datalist[[i+1]]<-impi(datalist[[i]])
      }
  }
  return(list(burnin=burninlist, data=datalist))
}

#' getProgramHelp
#' @export

getProgramHelp<- function ()
{
  paste("Program Arguments:\n",
        " Options:",
        "  -f (string)  Prior frequencies; comma separated pattern=frequency pairs, no spaces eg 00000=0.8,10000=0.05,01000=0.05",
        "  -s (integer) Program seed",
        "  -p Output Convergence Plot\n",
        " Required Arguments (6):",
        
        "  Path to Project Source File (./path/to/my/Rfile) (string)   required;  if set to the reserved string DEFAULT, use the sample package polio project",
        "  Number of Chains (integer)                                  required",
        "  Burn-in Reps (integer)                                      required;  if set to >0 do not monitor convergence during burn-in",
        "  Post-Burn-in Reps (integer)                                 required",
        "  Prior Sample Size Nprior (integer)                          required",
        "  Prefix for files generated                                  required",
        "\n", sep = "\n")
}


#' parsePriorFreqs
#' @param argstring string
#' @export
parsePriorFreqs<-function(argstring)
{
    v<-unlist(strsplit(argstring,","))
    l<-strsplit(v,"=")
    priorfreqs<-lapply(l, function(x) as.numeric(x[2]) )
    priorpatterns<-unlist(lapply(l, function(x) x[1]))
    if (length(priorfreqs) != length(priorpatterns)) stop("Error in parsePriorFreqs: prior frequencies must be specified as pattern1=freq1,pattern2=freq2,...etc")
    names(priorfreqs)<-priorpatterns
    priorfreqs
}

#' parseArgs
#' @export
parseArgs<-function()
{
	seed.default<-12112016
        n.nonopt<-6
                argv<-commandArgs(TRUE)
        nargv<-length(argv)
        if (argv[1] %in% c("-h","-?","-help","--h","--?","--help")) {
          cat(getProgramHelp())
          q()
        }
        if (nargv < n.nonopt) stop("Error in argument list\nUsage: program requires ", n.nonopt, " non-optional args:")
                nonopt<-argv[(nargv-n.nonopt +1) : nargv]
        if (sum(grepl("^-",nonopt,perl=T))>0) stop("Error in argument list\nUsage: program requires ",n.nonopt, " non-optional args:")
        	opt<-list()
	if (nargv > n.nonopt) {
          optval<-argv[1 : (nargv-n.nonopt)]
          loptval<-length(optval)
                    optpos<-grep("^-\\D", optval, perl=TRUE)
                    optval2<-optval
          optval2[optpos]<-sub("-", "", optval2[optpos])
                    for (i in optpos)
          {
                opt[[ optval2[i] ]]<- if (i == loptval || (i+1) %in% optpos) TRUE else
                        optval2[i+1]
          }
          if (length(nonopt) < n.nonopt) {
            cat("Usage: program requires ", n.nonopt, " non-optional args:")
            stop("Error in argument list\n")
          }
	}

        cat("\n***** Arguments *****\noptions:",paste(names(opt), opt),"\nrequired:",unlist(nonopt),"\n\n")

        return ( list(
                project=as.character(nonopt[[1]]),
                nchain=as.integer(nonopt[[2]]),
                burnin=as.integer(nonopt[[3]]),
                rep=as.integer(nonopt[[4]]),
                nprior=as.numeric(nonopt[[5]]),
                outprefix=as.character(nonopt[[6]]),
                fprior.list=if ('f' %in% names(opt)) parsePriorFreqs(opt[["f"]]) else NULL,
                seed=if ('s' %in% names(opt)) as.integer(opt[["s"]]) else seed.default,
                plot=if ('p' %in% names(opt)) TRUE else FALSE
        ))
}

#' loadProject
#' @param project.file project.file
#' @export
loadProject<-function(project.file)
{
	if (project.file %in% c("DEFAULT","polio","POLIO")) 
 	{
		loadProject("polioProject.R")
	} else {
 		if(!file.exists(project.file)) stop("Error in loadProject: Cannot find project file:",project.file,"\nmis-specified file name or path ?")
		source(project.file)
	}
}


#' nicePrintList
#' @param l length
#' @param col1 column1
#' @param col2 column2
#' @export
#' 
nicePrintList<-function(l,col1,col2)
{
	if (length(l)<1) return ()
        df<-data.frame(temp1=names(l),temp2=unlist(lapply(l,function(e) e[[1]])),row.names=NULL)
        names(df)<-c(col1,col2)
	print(df)
        cat("\n\n")
}

#' getChainSet
#' @param data Data frame or matrix of multivariate binary outcome with missing values to impute, NA denotes missing values
#' @param cols NULL or a character or integer vector giving the columns for imputation
#' @param nprior integer
#' @param fprior.list a list of prior frequencies. See "Details"
#' @param burnin number of burnin iterations, NULL or positive integer
#' @param rep number of iterations for imputation, default is 10
#' @return list of imputed data frames
#' @export
getChainSet<-function(data, cols=NULL, nchain=10, nprior, fprior.list, burnin, rep)
{
  sim.list<-list()

  if (missing(data)) ("Error in getChainSet: missing required 'data' argument\n")
  cat("\n***** Run Chains *****\n")
  for (d in 1:nchain){
    cat("now running: ", paste0(d, paste0("/", nchain)), "\n")	
    # impouta is a list of size 2 of lists:
    # the first element is a list of the burnin set, the second element is a list of size nreps
    # each rep is nind x noutcomes (T) 

    impouta<-milonga(data=data, cols=cols, nprior=nprior, fprior.list=fprior.list,
		burnin=burnin, rep=rep)

    # output the list of the rep set in the second element of the output list, size 2
    # burnin set in impouta[[1]] is not returned
    sim.list[[d]]<-impouta[[2]]
  }
  return(sim.list)
}

#' monitorConvergence
#' @param sim.list list of simulation ouput 
#' @param FUN function
#' @param ... other parameters
#' @export
monitorConvergence<-function(sim.list,FUN=monitorStats.default,...){
  FUN <- match.fun(FUN)
  nchain<-length(sim.list)
  nrep<-length(sim.list[[1]])
  stat<-NULL

  cat("\n***** Monitor Chains *****\n")
  r<-1
  cat("monitoring step: ", r, "\n")
  res.l<-lapply(sim.list, function ( this.chain,...) { FUN(this.chain[[r]], ...) },... )
  # number of variables to monitor
  nvarmon<-length(res.l[[1]])
  # res.arr is a 3D array
  res.arr<-array(0,c(nrep, nchain, nvarmon))
  # this gives an nchain x nvarmon matrix
  res.mat<-matrix(unlist(res.l), ncol=nvarmon, byrow=T)
  res.arr[r,,]<-res.mat
  for (r in 2:nrep)
  {
    cat("monitoring step: ", r, "\n")
    res.l<-lapply(sim.list, function ( this.chain,...) { FUN(this.chain[[r]], ...) },... )
    res.arr[r,,]<-matrix(unlist(res.l), ncol=nvarmon, byrow=T)
    stat<-rbind(stat, getCumulConvergeStat(res.arr[1:r,,,drop=FALSE]))
  }
  return(stat)
}

#' monitorStats.default
#' @param cmat matrix
monitorStats.default<-function(cmat)
{
    mean(ifelse(base::rowSums(cmat)>0, 1, 0) )
}

#' chainMIVars.default
#' @param cmat matrix
chainMIVars.default<-function(cmat)
{
    # prevalence of any positive  outcome
    prev<-mean(ifelse(base::rowSums(cmat)>0, 1, 0) )
    return( list(
       prev,
       prev*(1-prev)/nrow(cmat)
    )) 
}

#' getCumulConvergeStat
#' @param res.arr results object
#' @export
getCumulConvergeStat<-function(res.arr)
{
    r<-dim(res.arr)[1]
    nchain<-dim(res.arr)[2]
    nvar<-dim(res.arr)[3]
    out<-data.frame(nchain=nchain, nvar=nvar, step=r, stringsAsFactors=F)
    for (i in 1:nvar)
    {
      # mean of the monitored stat within chains
      psimean.d<-base::colMeans(res.arr[,,i,drop=FALSE],na.rm=TRUE)
      # var within each chain (d)
      psivar.d<-apply(res.arr[,,i,drop=FALSE], 2,var)
      # grand var  
      psivar<-var(psimean.d)

      #variance between
      B<-r*psivar*(1+nchain)/nchain
  
      #variance within
      W<-mean(psivar.d)

      # pooled var
      Vhat<-(W* (r-1) + B)/r
      sqrt_r<-sqrt(Vhat/W)
      suffix<-if (i==1) "" else as.character(i)
      nc<-ncol(out)
      out<-cbind(out,signif(B,4), signif(W,4), signif(Vhat,4), signif(sqrt_r,4))
      names(out)[(nc+1): ncol(out)] <- paste(c("varBetween","varWithin","varPooled","sqrtR"),suffix,sep="")
    }
    return(out)
}

#' getLastStep
#' @param sim.list list of simulation output
#' @export
#' 
getLastStep<-function(sim.list)
{
	nrep<-length(sim.list[[1]])
	lapply(sim.list, function (this.chain) { this.chain[[nrep]] })
}

#' runMITest
#' @param sim.list list of simulation output
#' @param alpha type I error
#' @param null.mean mean under NULL
#' @param FUN function
#' @param ... pass to other functions
#' @export
runMITest<- function (sim.last, alpha = 0.05, null.mean = 0, FUN = chainMIVars.default, ...)
{
  FUN <- match.fun(FUN)
  zcrit <- qnorm(alpha/2, lower = FALSE)
  nchain <- length(sim.last)
  mistat <- unlist(lapply(sim.last, FUN, ...))
  if (length(mistat) != 2 * nchain)
    stop("Error in runMITest: Need statistic + variance returned for MI inference for each complete data set\n")
  mistat.mat <- matrix(mistat, ncol = 2, byrow = T)
  D <- nrow(mistat.mat)
  mean.eta <- mean(mistat.mat[, 1])
  W <- mean(mistat.mat[, 2])
  B <- var(mistat.mat[, 1])
  Badj <- (D + 1) * B/D
  T <- W + Badj
  gamma <- Badj/T/df.k
  rel.var <- Badj/W
  df <- (D - 1)/gamma/gamma
  if (exists("df0") && !is.null(df0) && !is.na(df0)) {
    df.obs <- (1 - gamma) * df0 * (df0 + 1)/(df0 + 3)
    df <- 1/(1/df + 1/df.obs)
  }
  tstat <- (mean.eta - null.mean)/sqrt(T)
  p.value <- (1 - pt(abs(tstat), df)) * 2
  L95 <- mean.eta - qt(alpha/2, df, lower = F) * sqrt(T)
  U95 <- mean.eta + qt(alpha/2, df, lower = F) * sqrt(T)
  lambda <- (rel.var + 2/(df + 3))/(rel.var + 1)
  rel.eff <- 1/(1 + lambda/D)
  return(list(nchain = D, Vw = signif(W,4), Vb = signif(B,4),  Vt = signif(T,4),
              rel.var = signif(rel.var, 4),
              gamma = signif(gamma, 4),
              lambda = signif(lambda, 4),
              rel.eff = signif(rel.eff, 4),
              mean.eta = signif(mean.eta, 4),
              mean.l95 = signif(L95, 4), mean.u95 = signif(U95, 4),
              tstat = signif(tstat, 4),
              df = signif(df, 4),
              p.value = signif(p.value, 4)
  ))
}



#' getPatternPrior
#' @param pattern pattern matrix
#' @param fprior.list list of strings match to the pattern
#' @export
getPatternPrior<-function(pattern, fprior.list){
  K<-nrow(pattern)
  if(!is.null(fprior.list)) {
    prior.sum<-Reduce("+", fprior.list)
    if(prior.sum>1) stop("Error in getPatternPrior: sum of fprior.list frequencies must be <= 1")
  }
  pattern.names<-apply(pattern,1,function(p) paste0(p, collapse=""))
  pattern.prior<-NULL
  # if no prior frequencies, the default is 1/K
  if(is.null(fprior.list)){
    pattern.prior<-rep(1/K, K)
  } else{
    if(sum(!names(fprior.list)%in%pattern.names)>0)
      stop("Error in getPatternPrior: at least one pattern name in the fprior.list is incorrect")
    if(length(fprior.list)<K)
      ekr<-(1-prior.sum)/(K-length(fprior.list))
    for (name in pattern.names){
      if(name%in%names(fprior.list)){
        pattern.prior[name]<-fprior.list[[name]]
      } else {pattern.prior[name]<-ekr}
    }
  }
  return(pattern.prior)
}

#' plotConvergence
#' @param statdf statistics data frame
#' @param outprefix prefix for output
#' @param plot.var variables to plot
#' @export
plotConvergence<-function(statdf,outprefix,plot.var="sqrtR")
{
  library(ggplot2)
  grDevices::pdf(paste(outprefix,"pdf",sep="."),height=7, width=7,paper='special')
  ggplot2::theme_set(theme_bw())
  ggplot2::theme_update(axis.title = element_text(face="bold", size=14), axis.text = element_text(size=14))
  # find all plot.vars in the names of the df with suffixes (could be mulitple monitored vars)
  all.plot.var<-names(statdf)[names(statdf) %in% paste( plot.var, c("",1:100),sep="")]
  i<-0
  for (v in all.plot.var) {
    i<-i+1
    colvec<-as.character(cut(statdf[,v],c(-Inf,0.9999999999,1.2,5,+Inf),labels=c("blue","green","orange","red")))
    p<-ggplot2::ggplot(statdf, aes_string("step", v))+ggplot2::geom_point(col=colvec)+ggplot2::xlab("Step")
    if (plot.var =="sqrtR") {
       p<- p + if (i==1) ggplot2::ylab(expression(bold(sqrt("R")))) else if 
           (i==2) ggplot2::ylab(expression(bold(sqrt("R(2)")))) else if
           (i==3) ggplot2::ylab(expression(bold(sqrt("R(3)")))) else
	   ggplot2::ylab(expression(bold(sqrt("R(other)"))))
    }
    p<-p + ggplot2::geom_hline(yintercept = 1,linetype=3)
    print(p)
  }
  grDevices::dev.off()
}
