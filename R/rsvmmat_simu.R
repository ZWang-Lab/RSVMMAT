#' Simulation for RSVMMAT test
#'
#' This function use pre-defined parameters to make the simulation data for the RSVMMAT test (including type I and power test)
#'
#' @param n.sample Numeric, sample size, number of individuals
#' @param n.time Numeric, number of measurements for each individual
#' @param par List, the parameters for the phenotype traits, including covaraites and individual specific time dependent random effects
#' @param time_cov Logical variable, indicating whether time effect is included in phenotypic traits
#' @param gene.count Numeric, number of genes
#' @param intercept Logical variable, indicating whether intercept is used in phenotypic traits
#' @param causalpara List, the parameters for causal genes
#' @param power Logical variable, indicating whether include disease genes in the generated genes
#' @param phe.model String, the phenotype model, two optional values: 'logistic', 'liability'
#' @param oversampling String, the ascertainment scheme, three optional value: 'random', 'baseline', 'sum'
#' 
#' 
#' 
#' @return A list object is returned to be used as object for RSVMMAT test
#' @export

rsvmmat_simu<-function(n.sample=1000, n.time=5, par=list(), 
    time_cov = TRUE, gene.count = 10, intercept=TRUE, causalpara = list(), power = FALSE, phe.model = 'logistic', oversampling = "random")
{

    if(missing(par) || length(par)==0 )
    {
        par <- list(b0 = -1.8, b1 = 0.5, b2= 0.5, btime = 1, 
            sig.a = 0.8, sig.b =0.8, sig.e = 0.8, 
            rho=0.7);
        if(phe.model == "logistic"){
            par$b0 = -1.9
        } 
    }
  if(missing(causalpara) || length(causalpara)==0 ){
    causalpara=list(positive.ratio=1,causal.prop = 0.6,coef.causal = 0.12,tvarying.prop=1,invarcoef.causal=0.04)
  }
  

  if(power){

    #### Ascertainment sampling #####

    if(oversampling =="random")
    {
        snp.mat <- simu_snp_mat(n.sample*30, gene.count);
       snp.mat =  snp.mat[[1]]
      cat('* No ascertainment, random sampling:\n')
        phe <- simu.binary.phe.power(n.sample*30, n.time, par, causalpara,intercept, time_cov,  phe.model,snp.mat)
        index <- sample(1:(n.sample*30), n.sample)
        ID.select <- paste("id", index, sep = "")
        phe$y <- phe$y[index,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; snp.mat = snp.mat[index,]

      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")

    }else if(oversampling =="baseline")
    {
       snp.mat <- simu_snp_mat(n.sample*30, gene.count);
        snp.mat =  snp.mat[[1]]
       cat('* Ascertainment based on baseline:\n')
      phe <- simu.binary.phe.power(n.sample*30, n.time, par,causalpara, intercept, time_cov, phe.model,snp.mat);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      index.g1 <-which(phe$y[,1] ==1);
      index.g0 <- which(phe$y[,1]==0);
      index.select.0 <- index.g0[sample(n.sample*0.5)];
      index.select.1 <- index.g1[sample(n.sample*0.5)];
      index.select<-c(index.select.0, index.select.1)
      ID.select <- paste("id", index.select, sep = "")

      phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; snp.mat = snp.mat[index.select,]
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")

    }else if(oversampling == "sum"){
      snp.mat <- simu_snp_mat(n.sample*80, gene.count);
      snp.mat =  snp.mat[[1]]
      cat('* Ascertainment based on sum:\n')
      phe<- simu.binary.phe.power(n.sample*80, n.time, par, causalpara,intercept, time_cov,  phe.model,snp.mat);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      row.sum<-apply(phe$y, 1, sum)
      index.g2 <- which(row.sum==n.time);
      index.g0 <- which(row.sum==0);
      index.g1 <- which(row.sum<n.time& row.sum>0);

      index.select.0 <- index.g0[sample(n.sample*0.05)];
      index.select.1 <- index.g1[sample(n.sample*0.9)];
      index.select.2 <- index.g2[sample(n.sample*0.05)];
      index.select <- c(index.select.0, index.select.1, index.select.2);
      ID.select <- paste("id", index.select, sep = "")
     phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),];snp.mat = snp.mat[index.select,]
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")
    } else {
        print('Error! Choose oversampling methods from random, baseline and sum')
    }
  }else{
        snp.mat <- simu_snp_mat(n.sample, gene.count);

    #### Ascertainment sampling #####

    if(oversampling =="random")
    {
      cat('* No ascertainment, random sampling:\n')
        phe <- simu.binary.phe.null(n.sample*30, n.time, par, intercept, time_cov, phe.model);#changed 12/07
        index <- sample(1:(n.sample*30), n.sample)
        ID.select <- paste("id", index, sep = "")
        phe$y <- phe$y[index,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; 

      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")

    }else if(oversampling =="baseline")
    {
      cat('* Ascertainment based on baseline:\n')
      phe <- simu.binary.phe.null(n.sample*30, n.time, par, intercept, time_cov,  phe.model);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      index.g1 <-which(phe$y[,1] ==1);
      index.g0 <- which(phe$y[,1]==0);
      index.select.0 <- index.g0[sample(n.sample*0.5)];
      index.select.1 <- index.g1[sample(n.sample*0.5)];
      index.select<-c(index.select.0, index.select.1)
      ID.select <- paste("id", index.select, sep = "")

      phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),]; 
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")

    }else if(oversampling == "sum"){
      cat('* Ascertainment based on sum:\n')
      phe<- simu.binary.phe.null(n.sample*80, n.time, par, intercept, time_cov,  phe.model);
      cat('Disease Prevalence: ', apply(phe$y, 2, mean), "\n")
      row.sum<-apply(phe$y, 1, sum)
      index.g2 <- which(row.sum==n.time);
      index.g0 <- which(row.sum==0);
      index.g1 <- which(row.sum<n.time& row.sum>0);

      index.select.0 <- index.g0[sample(n.sample*0.05)];
      index.select.1 <- index.g1[sample(n.sample*0.9)];
      index.select.2 <- index.g2[sample(n.sample*0.05)];
      index.select <- c(index.select.0, index.select.1, index.select.2);
      ID.select <- paste("id", index.select, sep = "")
     phe$y <- phe$y[index.select,];phe$cov <- phe$cov[which(rownames(phe$cov)%in%ID.select),];
     # }
     cat("After Ascertainment: ", apply(phe$y, 2, mean), "\n")
    } else {
        print('Error! Choose oversampling methods from random, baseline and sum')
    }
  }  
  
    
 

    colnames(phe$y) <- paste("Y", 1:(NCOL(phe$y)), sep="") ;
    rownames(phe$y) <- paste("ID", 1:NROW(phe$y), sep="");

    colnames(phe$cov) <- paste("X", 1:(NCOL(phe$cov)), sep="") ;
    rownames(phe$cov) <- rep(paste("ID", 1:n.sample, sep=""), each = n.time);

    y.long <- c(t(phe$y));
    y.time <- cbind(rep(1:nrow(phe$y), each =n.time), rep((1:n.time)/n.time, nrow(phe$y)))
    # cov.long <- phe$cov;

    # rownames(cov.long) = NULL

    return(list(phe.wide = phe$y, phe.long = y.long, phe.time = y.time, phe.cov.long=phe$cov, snp.mat = snp.mat, phe.model = phe.model))
   
}

 # generate random effect
f.simu<-function( sample, n.time, par)
{
    ncol <- n.time;

    AR1 <- array(0,dim=c(ncol,ncol));
    for(i in 1:ncol)
    for(j in 1:ncol)
        AR1[i,j] <- par$rho^abs(i-j);

    sigma.b <- par$sig.b^2*AR1;
     
    r <- stats::rnorm( sample,  0, par$sig.a ) %*% array(1, dim=c(1,ncol))+
      mvtnorm::rmvnorm( sample,  rep(0, ncol), sigma.b ) ;
  
    return(r);
}



#generate null phenotype
simu.binary.phe.null<-function( n.sample, n.time, par, intercept, time_cov,  phe.model){

    cov.mat <- cbind( stats::rnorm(n.sample*n.time, 0, 1 ), rep(ifelse(stats::runif(n.sample)>0.5, 0, 1), each = n.time));

    if(intercept){
        mu <- f.simu(n.sample, n.time, par) + matrix(cbind(1, cov.mat )%*%c( par$b0, par$b1, par$b2 ), n.sample, n.time, byrow = TRUE)
    }else{
        mu <- f.simu(n.sample, n.time, par) + cov.mat %*% c(par$b1, par$b2 );
    }

    if(time_cov == TRUE){
        time.effect <- rep(1, n.sample)%*%(t(seq(1, n.time)*par$btime/n.time));
        mu = time.effect+mu;
    }


    #prepare output
    mu0=mu
    mu = c(t(mu));
    if(phe.model=='logistic')
    {
        print('logistic phenotypes')
        y <- matrix(stats::rbinom(n.sample*n.time, 1, inv.logit(mu)),  n.sample,n.time, byrow = T)
    }else if (phe.model=='liability'){#liability model
        print('liability phenotypes')
        #add sigma.e to the mu;
        # print(head(mu))
        mu <- mu + stats::rnorm( n.sample*n.time,  0, par$sig.e);
        # print(head(mu))
        y<- matrix(ifelse(mu>0, 1, 0), n.sample, n.time, byrow =T)
    }else {
        print('Wrong parameter, please choose from logistic, liability\n')
    }

    # cov.mat <- matrix(cov.mat, n.sample, n.time, byrow =T)

    rownames(cov.mat) <- rep(paste("id", 1:n.sample, sep=""), each = n.time);
    return(list(y=y, cov=cov.mat,mu=mu0, phe.model = phe.model));
}


inv.logit=function(x){
    return(exp(x)/(1+exp(x)))
}

logit=function(x){
    return(log(x/(1-x)))
}





simu_snp_mat<-function(  n.sample, mat.count)
{
  file.snp.hap1 <- system.file("extdata", "skat-test-1.hap.gz", package="RSVMMAT");
  file.snp.pos2 <- system.file("extdata", "skat-test-1.pos", package="RSVMMAT");
  
  snp.hap <- utils::read.table(file.snp.hap1, header=F);
  snp.pos <- utils::read.table(file.snp.pos2, header=T);
  
	snp.maxpos <- max(snp.pos$CHROM_POS);
	snp.minpos <- min(snp.pos$CHROM_POS);

	rare.cutoff <- 0.05;
      genelength <- 30*1000

	mat.ret <- list();
	n.mat <- 1;
	avoid_range <- c();
	while(n.mat<=mat.count)
	{
		snp.start <- as.integer(stats::runif(1, min=snp.minpos, max=snp.maxpos- genelength));
             snp.ends <- snp.start +  genelength;
		p.sel<-which( snp.pos$CHROM_POS>snp.start & snp.pos$CHROM_POS<snp.ends );


		snp.mat1 <- snp.hap[sample(nrow(snp.hap),n.sample,replace=T), p.sel + 2 ];
		snp.mat2 <- snp.hap[sample(nrow(snp.hap),n.sample,replace=T), p.sel + 2 ];
		snp.mat <- snp.mat1 + snp.mat2 -2;

		maf <- colMeans(snp.mat)/2;
		m.same <- which( maf==1 | maf==0 );
		if (length(m.same)>0)
			snp.mat <- snp.mat[, -m.same, drop=F ];

		#make sure of 2: Minor Allels, 0:Major Allels
		maf <- colMeans(snp.mat)/2;
		m.minor <- which( maf>0.5 );
		if (length(m.minor)>0)
			for(i in 1:length(m.minor))
				snp.mat[,m.minor[i]] <- 2 - snp.mat[,m.minor[i]];

		maf <- colMeans(snp.mat)/2;
		m.smallmaf <- which( maf < 0.001);
		if (length(m.smallmaf)>0)
			snp.mat <- snp.mat[, -m.smallmaf,drop=F ];

		if (dim(snp.mat)[2]==1)
			next;

		maf <- colMeans(snp.mat)/2;
		n.rare <- length(which(maf<rare.cutoff));
		if (n.rare<1)
			next;

		m.commonmaf <- which( maf > rare.cutoff)
            if( length(m.commonmaf)>0)
                 snp.mat <- snp.mat[, -m.commonmaf,drop=F ];


		rownames(snp.mat) <- paste("ID", 1:NROW(snp.mat), sep="");
		mat.ret[[n.mat]] <- snp.mat;

		n.mat <- n.mat + 1;
	}

	return(mat.ret);
}


simu_snp.ld<- function(n.sample, snp.count = 1000, ld = 0.5){
    ld.matrix = (1-ld)*diag(snp.count) + ld*rep(1, snp.count)%*%t(rep(1, snp.count))
    latent_cont1 = mvtnorm::rmvnorm(n.sample, rep(0, snp.count), ld.matrix)
    latent_cont2 = mvtnorm::rmvnorm(n.sample, rep(0, snp.count), ld.matrix)

    maf.rare = stats::runif(n = floor(snp.count*0.5), 0.01, 0.1 )
    maf.common = stats::runif(n= snp.count-floor(snp.count*0.5), 0.1, 0.5)
    maf = c(maf.rare, maf.common)

    snp.mat1 <-  snp.mat2  <- matrix(0, n.sample, snp.count)
    for(j in 1:snp.count){snp.mat1[,j] = ifelse(stats::pnorm(abs(latent_cont1[,j]), lower.tail = F)<maf[j]/2,1, 0)}
    for(j in 1:snp.count){snp.mat2[,j] = ifelse(stats::pnorm(abs(latent_cont2[,j]), lower.tail = F)<maf[j]/2,1, 0)}
    snp.mat <- snp.mat1 + snp.mat2 ;
    return(snp.mat)
}




simu.binary.phe.power  <- function(n.sample, n.time, par,causalpara, intercept, time_cov,  phe.model, snp.mat){
    y.null = simu.binary.phe.null(n.sample, n.time, par, intercept, time_cov,  phe.model)
    y0 = y.null$y; cov.mat = y.null$cov;  mu0 = y.null$mu
    maf <- colMeans(snp.mat)/2;
    nc = ncol(snp.mat);
    n.varyingcausal = round( causalpara$tvarying.prop*causalpara$causal.prop*nc)

    if(n.varyingcausal>0){
    gammat=logitcurve(n.time,causalpara$coef.causal)
    sign.causal = rep(0,n.varyingcausal)
        for(i in 1:n.varyingcausal){
            sign.causal[i] = ifelse(stats::runif(1)<causalpara$positive.ratio,1, -1)
            snp.effect <- sign.causal[i]*gammat*abs(log10(maf[i]))
            mu0 = mu0 +snp.mat[,i]%*% t(snp.effect)
        }
    }

    n.invariantcausal = round( (1-causalpara$tvarying.prop)*causalpara$causal.prop*nc)
    if(n.invariantcausal>0){
    gammat=rep( causalpara$invarcoef.causal,n.time)
    sign.causal = rep(0,n.invariantcausal)
        for(i in (n.varyingcausal+1):(n.varyingcausal+n.invariantcausal)){
            sign.causal[i] = ifelse(stats::runif(1)<causalpara$positive.ratio,1, -1)
            snp.effect <- sign.causal[i]*gammat*abs(log10(maf[i]))
            mu0 = mu0 +snp.mat[,i]%*% t(snp.effect)
        }
    }

    mu = c(t(mu0));
    if(phe.model=='logistic')
    {
        print('logistic phenotypes')
        y <- matrix(stats::rbinom(n.sample*n.time, 1, inv.logit(mu)),  n.sample,n.time, byrow = T)
    }else if (phe.model=='liability'){#liability model
        print('liability phenotypes')
        y<- matrix(ifelse(mu>0, 1, 0), n.sample, n.time, byrow =T)
    }else {
        print('Wrong parameter, please choose from logistic, liability\n')
    }

    return(list(y = y, cov= cov.mat, snp.mat = snp.mat))
}

logitcurve=function(ntime,coef.causal){
timeseq=(1:ntime)/(ntime)
ef=coef.causal/(1+coef.causal*exp(10-12*timeseq))
return(ef)
}



