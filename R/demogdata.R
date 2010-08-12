####################################################
### FUNCTIONS FOR HANDLING DEMOGRAPHIC RATE DATA ###
####################################################

## BASIC DEMOGDATA STRUCTURE:
# year - time series
# age - numeric vector of mid points of age groups
# rate - main rate (mortality, fertility or migration, etc)
#          matrices: female, male and total.
# pop - population matrices matching rates
# label (e.g., country, state, etc.) - character string
# type - mortality, fertility, migration, etc. character string.

demogdata <- function(data, pop, ages, years, type, label, name, lambda)
{
    p <- nrow(data)
    n <- ncol(data)
    if(nrow(pop) != p | ncol(pop) != n)
        stop("data and pop are of different size")
    if(length(ages) != p)
        stop("Number of ages doesn't match data")
    if(length(years) != n)
        stop("Number of years doesn't match data")

    types <- c("mortality","fertility","migration")
    idx <- pmatch(type,types)
    if(is.na(idx))
        warning("Unknown type")
    else
        type <- types[idx]

    if(missing(lambda))
    {
        if(type=="mortality")
            lambda <- 0
        else if(type=="fertility")
            lambda <- 0.4
        else
            lambda <- 1
    }

    obj <- list(year=years, age=ages, rate=list(as.matrix(data)), pop=list(as.matrix(pop)), type=type,
        label=label, lambda=lambda)
    dimnames(obj$rate[[1]]) <- dimnames(obj$pop[[1]]) <- list(ages,years)
    names(obj$rate) <- names(obj$pop) <- name
    return(structure(obj,class="demogdata"))
}

# Function to read demographic rates from file
# Assumed input format:
#  col 1: year
#  col 2: age
#  col 3+: rates for males, females and total (any order but labels required)
# Assume age cycles within year and all in order.
# Population file of same format.
# Output format:  object of class demogdata

read.demogdata <- function(file,popfile,type,label,max.mx=10,skip=2,popskip=skip,lambda)
{
    if(missing(lambda))
    {
        if(type=="mortality")
            lambda <- 0
        else if(type=="fertility")
            lambda <- 0.4
        else
            lambda <- 1
    }
    
    mfile <- !missing(file)
    mpopfile <- !missing(popfile)

    obj <- list(type=type,label=label,lambda=lambda)
    
    if(mfile)
    {
        tmp1 <- read.table(file,header=TRUE,na.strings=".",skip=skip)
        obj$year=sort(unique(tmp1[,1]))
        n <- length(obj$year)
        m <- length(unique(tmp1[,2]))
        obj$age <- tmp1[1:m,2]
        mnames <- names(tmp1)[-c(1,2)]
        n.mort <- length(mnames)
        obj$rate <- list()
        for(i in 1:n.mort)
        {
            obj$rate[[i]] = matrix(tmp1[,i+2], nrow=m, ncol=n)
            # Check bounds
            obj$rate[[i]][obj$rate[[i]] < 0] <- NA
            obj$rate[[i]][obj$rate[[i]] > max.mx] <- max.mx
            dimnames(obj$rate[[i]]) <- list(obj$age,obj$year)
        }
        names(obj$rate) = tolower(mnames)
    }

    if(mpopfile)
    {
        tmp2 <- read.table(popfile,header=TRUE,na.strings=".",skip=popskip)
        obj$year=sort(unique(tmp2[,1]))
        n <- length(obj$year)
        m <- length(unique(tmp2[,2]))
        obj$age <- tmp2[1:m,2]
        pnames <- names(tmp2)[-c(1,2)]
        if(mfile)
        {
            if(sum(pnames==mnames) != length(pnames))
                warning("Population names different from rates names")
            if(n!=ncol(obj$rate[[1]]) | m != nrow(obj$rate[[1]]))
                warning("Population matrices different size from rates matrices")
        }
        p.mort <- length(pnames)
        obj$pop <- list()
        for(i in 1:p.mort)
        {
            obj$pop[[i]] = matrix(tmp2[,i+2], nrow=m, ncol=n)
            # Check bounds
            obj$pop[[i]][obj$pop[[i]] < 0] <- NA
            dimnames(obj$pop[[i]]) <- list(obj$age,obj$year)
        }
        names(obj$pop) = tolower(pnames)
    }

    junk <- options(warn=-1)
    obj$age <- as.numeric(as.character(obj$age))
    options(warn=junk$warn)
    if(is.na(obj$age[m]))
        obj$age[m] <- 2*obj$age[m-1] - obj$age[m-2]

    return(structure(obj,class="demogdata"))
}

plot.demogdata <- function(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]), 
    datatype=ifelse(!is.null(x$rate),"rate","pop"),
    years=x$year, ages=x$age, max.age=max(x$age), transform=(x$type=="mortality"),
    plot.type= c("functions", "time", "depth", "density"), show.labels=TRUE, type="l", main=NULL, xlab, ylab,...)
{
    plot.type <- match.arg(plot.type)
    series <- tolower(series)
    ages <- ages[ages <= max.age]
    data <- extract.ages(extract.years(x,years),ages,FALSE)

    # Extract data matrix
    if(!is.element(datatype,names(data)))
        stop(paste("Data type",datatype,"not found"))
    tmp <- data[[pmatch(datatype,names(data))]]
    y <- get.series(tmp,series)

    # Transformation
	if(data$lambda > 1-1e-6)
		transform <- FALSE
    if(transform)
    {
        if(datatype=="rate")
            y <- BoxCox(y,data$lambda)
        else # Population
            y <- log(y)
    }

    # Choose appropriate y axis label
    if(missing(ylab))
    {
        if(datatype=="pop")
            ylab <- "Population"
        else if(data$type=="mortality")
            ylab <- "Death rate"
        else if(data$type=="fertility")
            ylab <- "Fertility rate"
        else if(data$type=="migration")
            ylab <- "Net migration"
        else if(data$type=="population")
            ylab <- "Population"
        else
            stop("This shouldn't happen!")
        if(transform)
        {
            if(data$lambda==0 | datatype=="pop")
                ylab <- paste("Log",tolower(ylab))
            else
                ylab <- paste("Transformed",tolower(ylab))
        }
    }

    # Choose appropriate axis title
    if(is.null(main))
    {
        if(data$type=="fertility" & series=="female")
            main <- data$label
        else
            main <- paste(data$label,": ",series,sep="")
        if(datatype=="pop")
            main <- paste(main,"population")
        else if(data$type=="mortality")
            main <- paste(main,"death rates")
        else if(data$type=="fertility")
            main <- paste(main,"fertility rates")
        else if(data$type=="migration")
            main <- paste(main,"net migration")
        if(length(years)==1)
            main <- paste(main,"  (",years,")",sep="")
        else
            main <- paste(main,"  (",min(years),"-",max(years),")",sep="")
    }

    # Produce plot
    if(length(data$age)==1)
    {
        main <- paste(main,"  Age:",data$age)
        if(missing(xlab))
            xlab <- "Year"
    }
    else if(missing(xlab))
        xlab <- "Age"

    plot(fts(data$age,y,s=years[1],f=1,yname="",xname=""),plot.type=plot.type,
        labels=switch(show.labels+1,NULL,rownames(y)),xlab=xlab,
        ylab=ylab,main=main,type=type,...)
}

lines.demogdata <- function(x, series=ifelse(!is.null(x$rate),names(x$rate)[1],names(x$pop)[1]), 
    datatype=ifelse(!is.null(x$rate),"rate",""),
    years=x$year, ages=x$age, max.age=max(x$age), transform=(x$type=="mortality"),
    plot.type= c("functions", "time", "depth", "density"), show.labels=TRUE,  ...)
{
    plot.type <- match.arg(plot.type)
    
    series <- tolower(series)
    ages <- ages[ages <= max.age]
    data <- extract.ages(extract.years(x,years),ages,FALSE)

    # Extract data matrix
    if(!is.element(datatype,names(data)))
        stop(paste("Data type",datatype,"not found"))
    tmp <- data[[pmatch(datatype,names(data))]]
    y <- get.series(tmp,series)

    # Transformation
    if(transform)
    {
        if(datatype=="rate")
            y <- BoxCox(y,data$lambda)
        else # Population
            y <- log(y)
    }

    # Set other arguments to appropriate values
    lines(fts(data$age,y,s=years[1],f=1),plot.type=plot.type, labels=switch(show.labels+1,NULL,rownames(y)),...)
}

points.demogdata <- function(...,pch=1)
{
    lines.demogdata(...,type="p",pch=pch)
}

print.demogdata <- function(x,...)
{
    Type <- x$type
    substr(Type,1,1) <- toupper(substr(Type,1,1))
    cat(paste(Type,"data for",x$label))
    cat("\n    Series: ")
    if(!is.null(x$rate))
        cat(names(x$rate))
    else
        cat(names(x$pop))
    cat(paste("\n    Years:",min(x$year),"-",max(x$year)))
    minx <- ifelse(min(x$age) < 0,"B",min(x$age))
    cat(paste("\n    Ages: ",minx,"-",max(x$age),"\n"))
}

summary.demogdata <- function(object, ...)
{
	print(object)
}
extract.years <- function(data,years)
{
    idx <- match(years,data$year)
    idx <- idx[!is.na(idx)]
	if(length(idx)==0)
		stop("No data available for those years")
    no.pop <- is.null(data$pop)
    no.rate <- is.null(data$rate)
    if(!no.rate)
    {
        nn <- length(data$rate)
        dname <- dimnames(data$rate[[1]])
        dname <- list(dname[[1]],dname[[2]][idx])
    }
    else if(!no.pop)
    {
        nn <- length(data$pop)
        dname <- dimnames(data$pop[[1]])
        dname <- list(dname[[1]],dname[[2]][idx])
    }
    else
        stop("No data!")
        
    for(j in 1:nn)
    {
        if(!no.rate)
        {
            data$rate[[j]] <- matrix(data$rate[[j]][,idx],ncol=length(idx))
            dimnames(data$rate[[j]]) <- dname
        }
        if(!no.pop)
        {
            data$pop[[j]] <- matrix(data$pop[[j]][,idx],ncol=length(idx))
            dimnames(data$pop[[j]]) <- dname
        }
        if(!is.null(data$obs.var))
        {
            if(length(dim(data$obs.var[[j]]))==2)
                data$obs.var[[j]] <- data$obs.var[[j]][,idx]
            else
                data$obs.var[[j]] <- data$obs.var[[j]][idx]
        }
    }
    data$year <- data$year[idx]
    return(data)
}


extract.ages <- function(data,ages,combine.upper=TRUE)
{
    if(combine.upper)
        data <- set.upperage(data,max(ages))

    idx <- match(ages,data$age)
    idx <- idx[!is.na(idx)]
    no.pop <- is.null(data$pop)
    no.rate <- is.null(data$rate)
    if(!no.rate)
    {
        nn <- length(data$rate)
        dname <- dimnames(data$rate[[1]])
        dname <- list(dname[[1]][idx],dname[[2]])
    }
    else if(!no.pop)
    {
        nn <- length(data$pop)
        dname <- dimnames(data$pop[[1]])
        dname <- list(dname[[1]][idx],dname[[2]])
    }
    else
        stop("No data!")

    for(j in 1:nn)
    {
        if(!no.rate)
        {
            data$rate[[j]] <- matrix(data$rate[[j]][idx,],nrow=length(idx))
            dimnames(data$rate[[j]]) <- dname
        }
        if(!no.pop)
        {
            data$pop[[j]] <- matrix(data$pop[[j]][idx,],nrow=length(idx))
            dimnames(data$pop[[j]]) <- dname
        }
        if(!is.null(data$obs.var))
        {
            if(length(dim(data$obs.var[[j]]))==2)
                data$obs.var[[j]] <- data$obs.var[[j]][idx,]
            else
                data$obs.var[[j]] <- data$obs.var[[j]][idx]
        }
    }
    data$age <- data$age[idx]

    return(data)
}


set.upperage <- function(data,max.age=100)
{
    if(max(data$age) < max.age)
        stop("max.age too large")
    else if(max(data$age) == max.age)
        return(data)

    if(is.null(data$pop))
        stop("This procedure needs the population data")
        
    no.rate <- is.null(data$rate)
    multiple.series <- is.list(data$pop)

    if(multiple.series)
    {
        nn <- length(data$pop)
        fred <- data
        for(j in 1:nn)
        {
            if(!no.rate)
                fred$rate <- data$rate[[j]]
            fred$pop <- data$pop[[j]]
            tmp <- set.upperage(fred,max.age)
            if(!no.rate)
                data$rate[[j]] <- tmp$rate
            data$pop[[j]] <- tmp$pop
        }
        data$age <- tmp$age
    }
    else
    {
        idx <- data$age >= max.age
        age <- data$age[!idx]
        rnames <- rownames(data$pop)[!idx]
        if(max(age) < max.age)
        {
            age <- c(age,max.age)
            rnames <- c(rnames,NA)
        }
        rnames[length(rnames)] <- paste(max.age,"+",sep="")
        upper.pop <- data$pop[idx,]
        pop = apply(matrix(upper.pop,nrow=sum(idx)),2,sum,na.rm=TRUE)
        data$pop <- rbind(matrix(data$pop[!idx,],ncol=ncol(data$pop)),pop)
        rownames(data$pop) <- rnames
        colnames(data$pop) <- data$year
        if(!no.rate)
        {
            upper.rate <- data$rate[idx,]
            actuals = apply(matrix(upper.rate*upper.pop,nrow=sum(idx)),2,sum,na.rm=TRUE)
            data$rate <- rbind(matrix(data$rate[!idx,],ncol=ncol(data$rate)),actuals/pop)
            rownames(data$rate) <- rnames
            colnames(data$rate) <- data$year
        }
        data$age <- age
        if(is.null(data$obs.var))
        {
            if(length(dim(data$obs.var))==2)
                data$obs.var <- data$obs.var[data$age <= max.age,]
            else
                data$obs.var <- data$obs.var[data$age <= max.age]
        }
    }
    return(data)
}


# Case insensitive version of is.element
is.el <- function(el,set)
{
    is.element(toupper(el),toupper(set))
}

get.series <- function(data,series)
{
    if(!is.el(series,names(data)))
        stop(paste("Series",series,"not found"))
    i <- match(toupper(series),toupper(names(data)))
    return(as.matrix(data[[i]]))
}

combine.demogdata <- function(obj1, obj2)
{
    if(!is.element("demogdata",class(obj1))  | !is.element("demogdata",class(obj2)))
        stop("Not demogdata objects")
    if(obj1$type != obj2$type)
        stop("Objects not of the same type")
    if(min(obj1$year) > min(obj2$year))
    {
        tmp <- obj2
        obj2 <- obj1
        obj1 <- tmp
    }
    if(max(obj1$year) > min(obj2$year))
        stop("Years overlap")

    pop <- (!is.null(obj1$pop) & !is.null(obj2$pop))
    idx <- match(obj1$age,obj2$age)
    idx <- idx[!is.na(idx)]
    m.obj <- list(label=obj1$label,age=obj1$age[idx],year=c(obj1$year,obj2$year),rate=list(),lambda = obj1$lambda)
    if(pop)
        m.obj$pop = list()
    if(!is.list(obj1$rate))
    {
        obj1$rate <- list(obj1$rate)
        if(pop)
            obj1$pop <- list(obj1$pop)
    }
    nn <- length(obj1$rate)
    k <- 0
    for (j in 1:nn)
    {
        i <- match(names(obj1$rate)[j], names(obj2$rate))
        if(is.na(i))
            i <- match(names(obj1$rate)[j], obj2$name)
        if(!is.na(i))
        {
            k <- k + 1
            m.obj$rate[[k]] <- cbind(obj1$rate[[j]][idx,],obj2$rate[[i]][idx,])
            names(m.obj$rate)[k] <- names(obj1$rate)[j]
            colnames(m.obj$rate[[k]]) <- m.obj$year
            rownames(m.obj$rate[[k]]) <- m.obj$age
            if(pop)
            {
                m.obj$pop[[k]] <- cbind(obj1$pop[[j]][idx,],obj2$pop[[i]][idx,])
                names(m.obj$pop)[k] <- names(obj1$pop)[j]
                colnames(m.obj$pop[[k]]) <- m.obj$year
                rownames(m.obj$pop[[k]]) <- m.obj$age
            }
        }
    }
    m.obj$type <- obj1$type
    return(structure(m.obj,class="demogdata"))
}

mean.demogdata <- function(x,series=names(x$rate)[1],transform=TRUE,na.rm=TRUE,...)
{
    mx <- get.series(x$rate,series)
    # Transformation
    if(transform)
        mx <- BoxCox(mx,x$lambda)

    mx[mx < -1e9] <- NA
    loc <- rowMeans(mx,na.rm=na.rm)
    return(list(x=x$age,y=loc))
}

median.demogdata <- function(x,series=names(x$rate)[1],
    transform=TRUE,method=c("hossjercroux","coordinate"),...)
{
    method = match.arg(method)
    mx <- get.series(x$rate,series)
    if(transform)
        mx <- BoxCox(mx,x$lambda)
    mx[mx < -1e9] <- NA
    loc <- L1median(t(mx),method=method)
    return(list(x=x$age,y=loc))
}

e0 <- function(data, series=NULL, years=data$year,
    type=c("period","cohort"), age=min(data$age), max.age=min(100,max(data$age)),
	PI=FALSE, nsim=500, ...)
{
    if(is.element("fmforecast",class(data)))
    {
		if(data$type != "mortality")
			stop("data not a mortality object")
        hdata <- list(year=data$model$year,age=data$model$age,
            type=data$type,label=data$model$label,lambda=data$lambda)
		if(min(data$model[[4]],na.rm=TRUE) > 0)
			hdata$rate <- list(data$model[[4]])
		else
		    hdata$rate <- list(InvBoxCox(data$model[[4]],data$lambda))
        names(hdata$rate) <- names(data$model)[4]
		if(!is.null(data$model$pop))
		{
		    hdata$pop = list(data$model$pop)
			names(hdata$pop) <- names(hdata$rate)
		}
        class(hdata) <- "demogdata"
        # Fix missing values. Why are they there?
        hdata$rate[[1]][is.na(hdata$rate[[1]])] <- 1-1e-5
        out <- structure(list(x=life.expectancy(hdata),
		    mean=life.expectancy(data,years=years,type=type,age=age,max.age=max.age),
            method="FDM model"),class="forecast")
		if(is.element("lca",class(data$model)))
			out$method = "LC model"
		else if(!is.null(data$product))
			out$method = "Coherent FDM model"
		if(PI) # Compute prediction intervals
		{
			if(is.null(data$product) & is.null(data$var) & is.null(data$kt.f))
				warning("Incomplete information. Possibly this is from a coherent\n  model and you need to pass the entire object.")
			else
			{
				sim <- simulate(data,nsim,...)
				e0sim <- matrix(NA,dim(sim)[2],dim(sim)[3])
				simdata <- data
				for(i in 1:dim(sim)[3])
				{
					simdata$rate[[1]] <- as.matrix(sim[,,i])
					e0sim[,i] <- life.expectancy(simdata)
				}
				if(is.element("lca",class(data$model)))
					out$level <- data$kt.f$level
				else
					out$level <- data$coeff[[1]]$level
				out$lower <- apply(e0sim,1,quantile,prob=0.5 - out$level/200)
				out$upper <- apply(e0sim,1,quantile,prob=0.5 + out$level/200)
				out$sim <- sim
			}
		}
		return(out)
    }
	else if(is.element("fmforecast2",class(data)))
	{
		if(data[[1]]$type != "mortality")
			stop("data not a mortality object")
		if(is.null(series))
			series <- names(data)[1]
		if(is.element("product",names(data))) # Assume coherent model
		{
			out <- e0(data[[series]],PI=FALSE)
			if(PI)
			{
				prodsim <- e0(data$product,nsim=nsim,PI=PI)
				ratiosim <- e0(data$ratio[[series]],nsim=nsim,PI=PI)
				#browser()
				sim <- prodsim$sim * ratiosim$sim
				e0sim <- matrix(NA,dim(sim)[2],dim(sim)[3])
				simdata <- data[[series]]
				for(i in 1:dim(sim)[3])
				{
					simdata$rate[[1]] <- as.matrix(sim[,,i])
					e0sim[,i] <- life.expectancy(simdata)
				}
				out$level <- data$product$coeff[[1]]$level
				out$lower <- apply(e0sim,1,quantile,prob=0.5 - out$level/200)
				out$upper <- apply(e0sim,1,quantile,prob=0.5 + out$level/200)
			}
		}
		else
			out <- e0(data[[series]],PI=PI,nsim=nsim)
		return(out)
	}		
    else
	{
	    if(!is.element("demogdata",class(data)))
			stop("data must be a demogdata object")
		if(data$type != "mortality")
			stop("data must be a mortality object")
		if(is.null(series))
			series <- names(data$rate)[1]
        return(life.expectancy(data,series=series,years=years,type=type,age=age,max.age=max.age))
	}
}



# Sex ratios

sex.ratio <- function(data) 
{
    if (class(data) == "demogdata") 
        rate.sr <- fts(x=data$age,y=data$rate$male/data$rate$female,start=min(data$year),f=1,xname="Age",yname="Year")
    else if(class(data) == "fmforecast2") 
        rate.sr <- fts(x=data$male$age,data$male$rate$male/data$female$rate$female,start=min(data$male$year),f=1,xname="Age",yname="Year")
    else
        stop("Unknown class of data")
    return(rate.sr)
}