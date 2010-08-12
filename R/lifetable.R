# Lifetable functions
# Produce lifetable from mortality rates
# Replicate the Excel spreadsheet life table
# This version only works with single ages.

lifetable <- function(data, series=names(data$rate)[1], years=data$year, ages=data$age,
    max.age=min(100,max(ages)), type=c("period","cohort"))
{
    if(!is.element("demogdata",class(data)))
        stop("data must be a demogdata object")
    if(data$type != "mortality")
        stop("data must be a mortality object")
    type <- match.arg(type)
    if(!is.el(series,names(data$rate)))
        stop(paste("Series",series,"not found"))
    if(is.na(sum(match(years,data$year))))
        stop("Years not present in data")
    if(max.age > max(ages))
        stop("max.age too large")
     
    sex <- series
    if(sex!="female" & sex!="male" & sex!="total")
    {
        if(is.element("model",names(data)))
            sex <- names(data$model)[4]
    }

    data <- extract.ages(data,ages,combine.upper=FALSE)
    if(max.age < max(ages))
        data <- extract.ages(data,min(ages):max.age,combine.upper=TRUE)
    if(type=="cohort")
        ltab = clifetab(data=data,series=series,years=years,sex=sex)
    else
        ltab = lifetab(data=data,series=series,years=years,sex=sex)

    ltab$label <- data$label

    idx <- match(ages,ltab$age)
    p <- length(idx)
    ltab$qx <- matrix(ltab$qx[idx,],nrow=p)
    ltab$mx <- matrix(ltab$mx[idx,],nrow=p)
    ltab$lx <- matrix(ltab$lx[idx,],nrow=p)
    ltab$dx <- matrix(ltab$dx[idx,],nrow=p)
    ltab$Lx <- matrix(ltab$Lx[idx,],nrow=p)
    ltab$Tx <- matrix(ltab$Tx[idx,],nrow=p)
    ltab$ex <- matrix(ltab$ex[idx,],nrow=p)

    rownames(ltab$qx) <- rownames(ltab$lx) <- rownames(ltab$dx) <- rownames(ltab$Lx) <- rownames(ltab$Tx) <- rownames(ltab$ex) <- ages

    if(type=="period")
    {
        ltab$rx <- matrix(ltab$rx[idx,],nrow=p)
        rownames(ltab$rx) <- ages-1
        rownames(ltab$rx)[ages==0] <- "B"
        colnames(ltab$qx) <- colnames(ltab$lx) <- colnames(ltab$dx) <- colnames(ltab$Lx) <- colnames(ltab$Tx) <- colnames(ltab$ex) <- colnames(ltab$rx) <- years[1:ncol(ltab$ex)]
    }
    else
    {
        ncol <- ncol(ltab$qx)
        colnames(ltab$qx) <- colnames(ltab$lx) <- colnames(ltab$dx) <- colnames(ltab$Lx) <- colnames(ltab$Tx) <- colnames(ltab$ex) <- paste("C",ages[1:ncol],sep="")
    }
    ltab$age = ages

    return(structure(ltab,class="lifetable"))
}


lifetab <- function(data,series,years,sex)
{
    subdata <- extract.years(data,years=years)
    mx <- get.series(subdata$rate,series)
    n <- length(years)
    p <- nrow(mx)
    rownames(mx) = subdata$age
    colnames(mx) = years
    qx = lx = dx = Lx = Tx = ex = rx = mx*NA
    rownames(rx) = subdata$age-1
    rownames(rx)[subdata$age==0] <- "B"
    startage = min(subdata$age)
    agegroup = subdata$age[4]-subdata$age[3]
    for(i in 1:n)
    {
        ltable <- lt(mx[,i],startage,agegroup,sex)
        p <- nrow(ltable)
        qx[1:p,i] = ltable$qx
        lx[1:p,i] = ltable$lx
        dx[1:p,i] = ltable$dx
        Lx[1:p,i] = ltable$Lx
        Tx[1:p,i] = ltable$Tx
        ex[1:p,i] = ltable$ex
        rx[1:p,i] = ltable$rx
    }
    
#    colnames(qx) <- colnames(lx) <- colnames(dx) <- colnames(Lx) <- colnames(Tx) <- colnames(ex) <- years
#    rownames(qx) <- rownames(lx) <- rownames(dx) <- rownames(Lx) <- rownames(Tx) <- rownames(ex) <- subdata$age
    return(list(age=subdata$age,year=years, mx=as.matrix(mx), qx=as.matrix(qx),
        lx=as.matrix(lx),dx=as.matrix(dx),Lx=as.matrix(Lx),Tx=as.matrix(Tx),ex=as.matrix(ex),rx=as.matrix(rx),
        series=series, type="period"))
}


# Cohort lifetable
#  Note number of years must be greater than number of ages
clifetab <- function(data,series,years,sex)
{
    idx <- match(years,data$year)
    idx <- idx[!is.na(idx)]
    if(length(idx)==0)
        stop("Year not available")
    years <- data$year[idx]
    mx <- get.series(data$rate,series)
    n <- ncol(mx)
    p <- xmx <- nrow(mx)
    if(n <= p)
    {
        p = n-1
        if(p<1)
            stop("Insufficient years")
        xmx <- nrow(mx)
#        mx <- as.matrix(mx[(1:p)+(xmx-p),]) # Is this right?
        mx <- as.matrix(mx[1:p,]) # Old method
    }
    idx <- idx[idx <= p]
    if(length(idx)==0)
        stop("Insufficient years")
    qx = dx = Tx = lx = Lx = ex = matrix(NA,p,p)
    ages <- data$age[(1:p)+(xmx-p)]
    rownames(lx) = rownames(Lx) = rownames(ex) = ages #is this right?
    colnames(lx) = colnames(Lx) = colnames(ex) = ages
    startage = min(ages)
    agegroup = ages[4]-ages[3]
    for (coh in 1:p)
    {
        cmx <- rep(0,p-coh+1)
        for (a in 1:(p-coh+1))
            cmx[a] <- mx[a+coh-1,a]
        ltable = lt(cmx, startage+coh-1,agegroup, sex=sex)
        lx[coh:p,coh] = ltable$lx[1:(p-coh+1)]
        Lx[coh:p,coh] = ltable$Lx[1:(p-coh+1)]
        ex[coh:p,coh] = ltable$ex[1:(p-coh+1)]
        qx[coh:p,coh] = ltable$qx[1:(p-coh+1)]
        dx[coh:p,coh] = ltable$dx[1:(p-coh+1)]
        Tx[coh:p,coh] = ltable$Tx[1:(p-coh+1)]
    }

    return(list(age=ages,year=years,mx=as.matrix(mx[,idx]),
        qx=as.matrix(qx[,idx]),lx=as.matrix(lx[,idx]),dx=as.matrix(dx[,idx]),
            Lx=as.matrix(Lx[,idx]),Tx=as.matrix(Tx[,idx]),ex=as.matrix(ex[,idx]),series=series,
        type="cohort"))
}

# Derive lifetable values from single year mortality rates
# mx is a vector consisting of one column from a mortality matrix

#oldlt <- function(mx,startage=0)
#{
#    # Omit missing ages
#    firstmiss <- (1:length(mx))[is.na(mx)][1]
#    if(!is.na(firstmiss))
#        mx <- mx[1:(firstmiss-1)]
#    nn <- length(mx)
#    if(startage==0)
#        a0 <- 0.07 + (1.7*mx[1]) # per Keyfitz
#    else
#        a0 <- 0.5
#    ax <- c(a0, rep(0.5, nn-1))
#    qx <- mx/(1 + ((1-ax) * mx)) ##Chiang
#    qx[nn]=1
#    qx[qx > 1] <- 1
#    lx <- c(1,cumprod(1-qx[1:(nn-1)]))
#    dx <- -diff(c(lx,0))
#    Lx <- lx - dx  + dx*ax
#    Lx[nn] <- lx[nn]/mx[nn]
#    Tx <- rep(0,nn)
#    Tx[nn] <- Lx[nn]
#    Tx <- rev(cumsum(rev(Lx)))
#    ex <- Tx/lx
#    rx <- c(Lx[1]/lx[1],Lx[2:(nn-1)]/Lx[1:(nn-2)],Tx[nn]/Tx[nn-1])
#    
#    result <- data.frame(mx = mx, qx = qx, lx = lx, dx = dx,
#                Lx = Lx, Tx = Tx, ex = ex, rx=rx)
#    return(result)
#}
#

lt <- function (mx, startage=0, agegroup=5, sex)
{
    # Omit missing ages
    if(is.na(mx[1]))
        mx[1] <- 0
    firstmiss <- (1:length(mx))[is.na(mx)][1]
    if (!is.na(firstmiss)) 
        mx <- mx[1:(firstmiss - 1)]
    nn <- length(mx)

    # Compute width of each age group
    if(agegroup == 1)
        nx <- c(rep(1,nn-1),Inf)
    else if(agegroup == 5) # First age group 0, then 1-4, then 5-year groups.
        nx <- c(1,4,rep(5,nn-3),Inf)     
    else
        stop("agegroup must be either 1 or 5")

    if (startage == 0) # for single year data and the first age (0) in 5-year data
    {
        if (sex == "female")
        {
            if (mx[1] < 0.107) 
                a0 <- 0.053 + 2.800 * mx[1] 
            else 
                a0 <- .35 
        }
        else if (sex == "male")  
        {
            if (mx[1] < 0.107) 
                a0 <- 0.045 + 2.684 * mx[1] 
            else 
                a0 <- .33 
        }
        else #if (sex == "total") 
        {
            if (mx[1] < 0.107) 
                a0 <- 0.049 + 2.742 * mx[1] 
            else 
                a0 <- .34 
        }
            
    }
    else if (startage > 0) 
        a0 <- 0.5
    if (agegroup == 1) 
        ax <- c(a0, rep(0.5, nn - 2),Inf) ## last ax not used 
    else if(agegroup == 5) 
    {
        if (sex == "female") 
        {
            if (mx[1] < 0.107) 
                a1 <- 1.522 - 1.518 * mx[1] 
            else 
                a1 <- 1.361 
        }
        else if (sex == "male") 
        {
            if (mx[1] < 0.107) 
                a1 <- 1.651 - 2.816 * mx[1] 
            else 
                a1 <- 1.352 
        }
        else #if (sex == "total") 
        {
            if (mx[1] < 0.107) 
                a1 <- 1.5865 - 2.167 * mx[1] 
            else 
                a1 <- 1.3565 
        }
        ax <- c(a0,a1,rep(2.6, nn-3),Inf)               
    }    ### ax=2.5 known to be too low esp at low levels of mortality
    if (startage > 0)  
    { 
        a0  <- 2.6
        if (agegroup == 5) 
        {
            ax <- c(a0,rep(2.6, nn-2),Inf)
            nx <- c(rep(5,nn))  
        }    ## last ax not used      
    }        
    
    qx <- nx * mx /(1 + (nx - ax)* mx)

    age <- startage + cumsum(nx) - 1
    
    if (max(age) >= 75)
    {
        idx <- age>=75
        ax[idx] <- (1/mx + nx - nx/(1-exp(-nx*mx)))[idx]
        qx[idx] <- 1 - exp(-nx * mx)[idx]
    }
    
    
    qx[qx > 1] <- 1  # this problem  avoided by recalc at 75+ ...should not be needed
    qx[nn] <- 1 
    lx <- c(1, cumprod(1 - qx[1:(nn - 1)]))
    dx <- -diff(c(lx, 0))
    Lx <- nx * lx - dx * (nx - ax)   # nx used here
    Lx[nn] <- lx[nn]/mx[nn]
    Tx <- rep(0, nn)
    Tx[nn] <- Lx[nn]
    Tx <- rev(cumsum(rev(Lx)))
    ex <- Tx/lx
    rx <- c(Lx[1]/lx[1], Lx[2:(nn-1)]/Lx[1:(nn-2)], Tx[nn]/Tx[nn-1])
    if( agegroup == 5) 
        rx <- c(0, (Lx[1]+Lx[2])/5*lx[1], Lx[3]/(Lx[1]+Lx[2]), 
            Lx[4:(nn-1)]/Lx[3:(nn-2)], Tx[nn]/Tx[nn-1]) # note rx[1]= 0 (not used) to maintain length = nn 
    result <- data.frame(ax = ax, mx = mx, qx = qx, lx = lx, dx = dx,
        Lx = Lx, Tx = Tx, ex = ex, rx = rx, nx =nx)
    return(result)
}


# Compute expected age from single year mortality rates
get.e0 <- function(x,agegroup,sex)
{
    lt(x, 0, agegroup, sex)$ex[1]
}

# Compute expected ages for multiple years from period lifetable
life.expectancy <- function(data,series=names(data$rate)[1],years=data$year,
    type=c("period","cohort"), age=min(data$age), max.age=min(100,max(data$age)))
{
    type <- match.arg(type)
    if(!is.el(series,names(data$rate)))
        stop(paste("Series",series,"not found"))
    if(age > max.age | age > max(data$age))
        stop("age is greater than maximum age")
    else if(age < min(data$age))
        stop("age is less than minimum age")
        
    data.lt <- lifetable(data,series,years,type=type,max.age=max.age)$ex
    idx <- match(age,rownames(data.lt))
    if(sum(is.na(data.lt[idx,]))>0 | max(data.lt[idx,]) > 1e9)
        warning("Some missing or infinite values in the life table calculation.\n  These can probably be avoided by setting max.age to a lower value.")

    return(ts(data.lt[idx,],s=years[1],f=1))
}

plot.lifetable <- function(x,years=x$year,main,xlab="Age",ylab="Expected number of years left",...)
{
    # Extract years
    idx <- match(years,x$year)
    idx <- idx[!is.na(idx)]
    idx <- idx[idx <= ncol(x$ex)]
    if(length(idx)==0)
        stop("Year not available")
    years <- x$year[idx]
    ny <- length(years)

    if(missing(main))
    {
        main <- paste("Life expectancy:",x$label,x$series)
        if(ny>1)
            main <- paste(main," (",min(years),"-",max(years),")",sep="")
        else
            main <- paste(main," (",years,")",sep="")
    }

    plot(fts(x$age,x$ex[,idx],s=years[1],f=1),main=main,ylab=ylab,xlab=xlab,...)
}

lines.lifetable <- function(x,years=x$year,...)
{
    # Extract years
    idx <- match(years,x$year)
    idx <- idx[!is.na(idx)]
    idx <- idx[idx <= ncol(x$ex)]
    if(length(idx)==0)
        stop("Year not available")
    years <- x$year[idx]
    ny <- length(years)

    lines(fts(x$age,x$ex[,idx],s=x$year[1],f=1),...)
}

print.lifetable <- function(x,years=x$year,ages=x$age,digits=4,...)
{
    # Extract years
    idx <- match(years,x$year)
    idx <- idx[!is.na(idx)]
    idx <- idx[idx <= ncol(x$ex)]
    if(length(idx)==0)
        stop("Year not available")
    years <- x$year[idx]
    ny <- length(years)
    outlist <- vector(length=ny,mode="list")

    # Extract ages
    idx3 <- match(ages,x$age)
    idx3 <- idx3[!is.na(idx)]
    if(length(idx3)==0)
        stop("Age not available")
    ages <- x$age[idx3]
    cohort <- colnames(x$ex)

    # Construct output
    for(i in 1:ny)
    {
        j <- idx[i]
        idx2 <- !is.na(x$ex[,j])
        idx2 <- idx2
        if(sum(idx2)>0)
        {
            outlist[[i]] <- data.frame(x$mx[,j],x$qx[,j],x$lx[,j],x$dx[,j],x$Lx[,j],x$Tx[,j],x$ex[,j])[idx2,]
            rownames(outlist[[i]]) <- rownames(x$ex)[idx2]
            colnames(outlist[[i]]) <- c("mx","qx","lx","dx","Lx","Tx","ex")
            idx4 <- match(ages,rownames(outlist[[i]]))
            idx4 <- idx4[!is.na(idx4)]
            if(length(idx4)==0)
                stop("Insufficient data")
            outlist[[i]] <- outlist[[i]][idx4,]
        }
    }
    names(outlist) = years
    if(x$type=="period")
        cat("Period ")
    else if(x$type=="cohort")
        cat("Cohort ")
    else
        stop("Unknown lifetable type")
    cat(paste("lifetable for",x$label,":",x$series,"\n\n"))
    for(i in 1:ny)
    {
        if(!is.null(outlist[[i]]))
        {
            if(x$type=="period")
                cat(paste("Year:",years[i],"\n"))
            else
                cat(paste("Cohort:",cohort[i],"\n"))
            print(round(outlist[[i]],digits=digits))
            cat("\n")
        }
    }
    invisible(outlist)
}
