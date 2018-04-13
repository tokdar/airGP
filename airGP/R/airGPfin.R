airGP <- function(x.train, y.train, nsamp = 100, thin = 10,
                  sparsity = list(ncomp = NULL, dmax = 4, Esize = 1, pactive = c(1,1),
                                  shrink = c("shape", "rate"), decay.type = c("poly", "geo")),
                  lowrank = list(max.rank = nrow(x.train), tol = 1e-6, lpenalty = NULL), 
                      hyper = list(ERsq = 0.95, lam = c(9.9, 0.1), rho = c(1, 5.0), sig = c(0,0),
                                   prox = c(0.99, 0.93, 0.75, 0.45, 0.15), nrho = 10),
                  adapt = list(var.p = NULL, init.meth = c("corr", "lasso", "gp1"), 
                               taper.imp = 0.95, simMat = NULL, simFn = cor), 
                  state = NULL, verbose = TRUE, repeat.call = FALSE){

    cl <- match.call()
    if(repeat.call) print(cl)
    #==================== DATA PREPROCESSING ==========================
    #-------------------- data dimension and predictor names ----------
    n <- length(y.train)
    xnames <- dimnames(x.train)[[2]]
    if(!is.matrix(x.train)) x.train <- as.matrix(x.train, nrow = n)
    p <- ncol(x.train)
    if(is.null(xnames)) xnames <- paste0("x", 1:p)
    #-------------------- translate each predictor to [0,1] -----------
    mx <- apply(x.train, 2, min)
    sx <- apply(x.train, 2, function(z) diff(range(z)))
    x <- scale(x.train, center = mx, scale = sx)
    dimnames(x) <- list(1:n, xnames)
    #-------------------- standardize response ------------------------
    my <- mean(y.train)
    sy <- sd(y.train)
    obs <- (y.train - my) / sy
    
    #===================== MODEL SETTING ==============================
    #--------------------- GP HYPERPARAMETERS -------------------------
    if(is.null(hyper$ERsq)) hyper$ERsq <- 0.95
    if(is.null(hyper$lam)) hyper$lam <- c(9.9, 0.1)
    if(is.null(hyper$rho)) hyper$rho <- c(1, 5.0)
    if(is.null(hyper$sig)) hyper$sig <- c(0, 0)
    if(is.null(hyper$nrho)) hyper$nrho <- 10
    
    lam.list <- lengthscale.fixation(hyper)
    lamsq <- lam.list$lamsq
    lplam <- lam.list$lplam
    
    
    #--------------------- SPARSITY CONTROL ---------------------------
    if(is.null(sparsity$ncomp)) sparsity$ncomp <- min(max(5, 2*ceiling(n/log(p, 10))), ceiling(sqrt(max(n,p))), p, 50)
    ncomp <- sparsity$ncomp
    dmax <- sparsity$dmax; if(is.null(dmax)) dmax <- 4; dmax <- min(dmax, p)
    Esize <- sparsity$Esize; if(is.null(Esize)) Esize <- 1
    pactive <- sparsity$pactive; shrink <- sparsity$shrink[1]; decay.type <- sparsity$decay.type[1]
    if(is.null(shrink)) shrink <- "rate"
    if(is.null(decay.type)) decay.type <- "poly"
    if(is.null(pactive)) pactive <- c(1,1)

    pactive.par <- beta.pars(pactive, c(0.5, 0.9999))$par
    lpmodelsize <- dpois(1:dmax, Esize, log = TRUE)
    lpmodelsize <- lpmodelsize - logsum(lpmodelsize) - lchoose(p, 1:dmax)

    rhoA.factor <- rep(1, ncomp)
    rhoB.factor <- rep(1, ncomp)
    if(shrink == "shape") {
        if(decay.type == "poly")
            rhoA.factor <- 1/(1:ncomp)^2
        else if (decay.type == "geo")
            rhoA.factor <- (0.5)^(1:ncomp)
        else
            rhoA.factor <- rep(1/ncomp, ncomp)
    } else if (shrink == "rate") {
        if(decay.type == "poly")
            rhoB.factor <- (1:ncomp)^2
        else if(decay.type == "geo")
            rhoB.factor <- (2.0)^(1:ncomp)
        else
            rhoB.factor <- rep(ncomp, ncomp)
    }
    #rhob.mean <- hyper$rho * (1 - hyper$ERsq) * sum(rhoA.factor / rhoB.factor) / hyper$ERsq
    rhob.mean <- hyper$rho[1] * (1 - hyper$ERsq) / hyper$ERsq
    rho.par <- c(hyper$rho, rhob.mean)
    
    #--------------------- LOW RANK APPROX CONTROL --------------------
    if(is.null(lowrank$max.rank)) lowrank$max.rank <- n
    if(lowrank$max.rank > n) lowrank$max.rank <- n
    if(is.null(lowrank$tol)) lowrank$tol <- 1e-6
    if(is.null(lowrank$lpenalty)) lowrank$lpenalty <- -dt(scale(x), df = 10, log = TRUE)
    
    
    ##-------------------- MCMC ADAPTATION CONTROL --------------------
    if(is.null(adapt$init.meth)) adapt$init.meth <- "corr"
    if(is.null(adapt$taper.imp)) adapt$taper.imp <- 0.95
    if(is.null(adapt$simFn)) adapt$simFn <- cor
    
    taper.imp <- max(min(1, adapt$taper.imp), 0)
    xx.sim <- adapt$simMat
    if(is.null(xx.sim)) xx.sim <- abs(adapt$simFn(x.train))
    taper.cut <- as.numeric(quantile(xx.sim[upper.tri(xx.sim)], taper.imp))
    
    nhbr <- lapply(1:p, function(i) which(xx.sim[,i] > taper.cut))
    nhbr.len <- sapply(nhbr, length)
    nhbr.val <- unlist(lapply(1:p, function(i) return(fnB(xx.sim[nhbr[[i]],i], taper.cut))))
    nhbr.ix <- unlist(nhbr) - 1 ## for c coding, index starts at 0
    
    var.p <- adapt$var.p
    if(is.null(var.p)){
        if(adapt$init.meth[1] == "corr"){
            xy.cor <- abs(apply(x.train, 2, function(z) cor(z, y.train)))
            var.p <- exp(xy.cor)
        } else if(adapt$init.meth == "lasso"){
            o.lasso <- lapply(1:4, function(j) cv.glmnet(x.train, y.train))
            beta.lasso <- rowMeans(sapply(o.lasso, function(o) return(coef(o, s = "lambda.min")[-1])))
            imp.lasso <- abs(beta.lasso)
            var.p <- rep(0, p)
            for(j in 1:p) var.p[nhbr[[j]]] <- var.p[nhbr[[j]]] + imp.lasso[j] * fnB(xx.sim[nhbr[[j]],j], taper.cut)
             var.p <- 1 + var.p / (1e-6/p + max(var.p))
        } else if(adapt$init.meth == "gp1"){
            fit.1v <- airGP1v(x.train, y.train, lowrank = list(max.rank = 21))
            imp.1v <- exp(fit.1v$ls1 - max(fit.1v$ls1))
            var.p <- rep(0, p)
            for(j in 1:p) var.p[nhbr[[j]]] <- var.p[nhbr[[j]]] + imp.1v[j] * fnB(xx.sim[nhbr[[j]],j], taper.cut)
            var.p <- 1 + 2 * var.p / (max(var.p) + 1e-6/p)
        } else {
            stop("adapt$init.meth must match 'corr', 'lasso' or 'gp1'")
        }
    }
    
    
    ## setting Markov chain sampler controls
    ##cat("Starting MCMC with summary(var.p) =", round(summary(var.p), 2), "\n")
    p.ref <- 0.0 * (0:p > 0) * (0:p <= dmax)
    cut.point <- min(dmax, p - 1)
    p.swap <- 2 * (1 - p.ref) * dpois(0:p, 3) * (0:p > 0) * (0:p <= cut.point)
    cut.point <- min(dmax, p)
    p.add <- (1 - p.ref) * dexp(0:p, 1) * (0:p < cut.point)
    p.rem <- 1 - p.add - p.swap - p.ref
    
    p.move <- rbind(p.add, p.rem, p.swap, p.ref)
    ##print(p.move[,1:min(10,p+1)])
    
    ## set initial state if asked for
    initiate.state <- is.null(state)
    if(initiate.state){
        p.active <- pactive.par[1] / sum(pactive.par)
        sigma <- 1
        rho.b <- rhob.mean
        comp.size <- rep(0, ncomp)
        comp.ix <- rep("", ncomp)
        gp.par <- rep(0, 2*ncomp)
        fcomp <- rep(0, ncomp*n)
    } else {
        p.active <- state$p.active; if(is.null(p.active)) stop("'state' missing 'pactive'")
        sigma <- state$sigma; if(is.null(sigma)) stop("'state' missing 'sigma'")
        rho.b <- state$rho.b; if(is.null(rho.b)) stop("'state' missing 'rho.b'")
        comp.size <- state$comp.size; if(is.null(comp.size)) stop("'state' missing 'comp.size'")
        comp.ix <- state$comp.ix; if(is.null(comp.ix)) stop("'state' missing 'ix'")
        gp.par <- state$gp.par; if(is.null(gp.par)) stop("'state' missing 'gp.par'")
        fcomp <- state$fcomp; if(is.null(fcomp)) stop("'state' missing 'fcomp'")
    }
    
    ##cat("Call C code\n")
    st <- system.time(
    out.c <- .C("airGP_tog", data = as.double(c(obs, c(x))), state = as.integer(initiate.state), state.pars = as.double(c(p.active, sigma, rho.b, gp.par)), fcomp = as.double(fcomp), comp.size = as.integer(comp.size), comp.ix = as.character(comp.ix), dim = as.integer(c(n, p, ncomp, lowrank$max.rank, length(lamsq), hyper$nrho, dmax, nsamp, thin)), nhbr.ix = as.integer(nhbr.ix), nhbr.len = as.integer(nhbr.len), nhbr.val = as.double(nhbr.val), gp.grid = as.double(c(lamsq, lplam)), hpar = as.double(c(hyper$sig, pactive.par, rho.par, hyper$lam)), comp.adj = as.double(c(rhoA.factor, rhoB.factor)), lpmodelsize = as.double(lpmodelsize), pmove = as.double(p.move), move.diagnostics = double(8), lpen = as.double(lowrank$lpenalty), varp = as.double(var.p), ix.store = character(nsamp * ncomp), par.store = double(2*ncomp*nsamp), sigma.store = double(nsamp), ftot.store = double(nsamp*n), pactive.store = double(nsamp), rhob.store = double(nsamp), tolchol = as.double(lowrank$tol), verbose = as.integer(verbose))
    )
    
    last.state <- list(p.active = out.c$state.pars[1], sigma = out.c$state.pars[2], rho.b = out.c$state.pars[3], comp.size = out.c$comp.size, comp.ix = out.c$comp.ix, gp.par = out.c$state.pars[3 + 1:(2*ncomp)], fcomp = matrix(out.c$fcomp, ncol = ncomp))
    
    out <- list(ix.store = t(matrix(out.c$ix.store, nrow = ncomp)), pars.store = t(matrix(out.c$par.store, nrow = 2 * ncomp)), sigma.store = out.c$sigma.store, ftot.store = matrix(out.c$ftot.store, nrow = n), pactive.store = out.c$pactive.store, rhob.store = out.c$rhob.store, state = last.state, nprop = out.c$move.diagnostics[1:4], nacpt = out.c$move.diagnostics[5:8], var.p = out.c$varp, my = my, sy = sy, mx = mx, sx = sx, x = x, obs = obs, hyper = hyper, lamsq.grid = lamsq, runtime = st[3], lowrank = lowrank, details = list(nhbr.ix = out.c$nhbr.ix, nhbr.len = out.c$nhbr.len, nhbr.val = out.c$nhbr.val, gp.grid = out.c$gp.grid, hpar = out.c$hpar, pmove = out.c$pmove, move.diagnostics = out.c$move.diagnostics, dim = out.c$dim, verbose = out.c$verbose)
    )
    class(out) <- "airGP"
    return(out)
}

airGP1v <- function(x.train, y.train, 
                    lowrank = list(max.rank = nrow(x.train), tol = 1e-6, lpenalty = NULL), 
                    hyper = list(ERsq = 0.95, lam = c(9.9, 0.1), rho = c(1,1), sig = c(0,0), delta = 0.1, 
                                 prox = c(0.99, 0.93, 0.75, 0.45, 0.15))){
    
    xnames <- dimnames(x.train)[[2]]
    
    n <- length(y.train)
    x.train <- matrix(x.train, nrow = n)
    p <- ncol(x.train)
    
    if(is.null(hyper$ERsq)) hyper$ERsq <- 0.95
    if(is.null(hyper$lam)) hyper$lam <- c(9.9, 0.1)
    if(is.null(hyper$rho)) hyper$rho <- c(1, 1)
    if(is.null(hyper$sig)) hyper$sig <- c(0, 0)
    if(is.null(hyper$delta)) hyper$delta <- 0.1
    
    mx <- apply(x.train, 2, min)
    sx <- apply(x.train, 2, function(z) diff(range(z)))
    x <- scale(x.train, center = mx, scale = sx)
    dimnames(x)[[2]] <- xnames
    
    if(is.null(lowrank$max.rank)) lowrank$max.rank <- n
    if(lowrank$max.rank > n) lowrank$max.rank <- n
    if(is.null(lowrank$tol)) lowrank$tol <- 1e-6
    if(is.null(lowrank$lpenalty)) lowrank$lpenalty <- -dt(scale(x), df = 10, log = TRUE)
    
    my <- mean(y.train)
    sy <- sd(y.train)
    obs <- as.numeric(scale(y.train))
    
    lam.list <- lengthscale.fixation(hyper)
    lamsq <- lam.list$lamsq
    lplam <- lam.list$lplam
    
    rho.list <- magnitude.fixation(hyper)
    rhosq <- rho.list$rhosq
    lprho <- rho.list$lprho
    
    lpmodelsize <- -lchoose(p, 1)
    
    st <- system.time(
    out.c <- .C("airGP_one", 
                xvar = as.double(x), 
                yvar = as.double(obs), 
                dim = as.integer(c(n, p, lowrank$max.rank, length(lamsq), length(rhosq))),
                lamsqR = as.double(lamsq), 
                rhosqR = as.double(rhosq), 
                hpar = as.double(hyper$sig), 
                lpmodelsize = as.double(lpmodelsize), 
                lplamR = as.double(lplam), 
                lprhoR = as.double(lprho), 
                lpen = as.double(lowrank$lpenalty), 
                ls1 = double(p), 
                ls0 = double(1), 
                tolchol = as.double(lowrank$tol))
    )
    return(list(ls1 = out.c$ls1, ls0 = out.c$ls0, runtime = st[3]))
}

predict.airGP <- function(object, x.test, burn = 0.1, nmc = 200, ...){
    
    ix.store <- object$ix.store
    pars.store <- object$pars.store
    
    nsamp <- nrow(pars.store)
    ncomp <- ncol(ix.store)
    
    n.postburn <- floor(nsamp * (1 - burn))
    nmc <- min(nmc, n.postburn)
    thin <- floor(n.postburn / nmc)
    ss <- nsamp - thin * ((nmc - 1):0)
    nss <- length(ss)
    
    x <- object$x
    obs <- object$obs
    
    n <- nrow(x)
    p <- ncol(x)
    xnew <- scale(x.test, center = object$mx, scale = object$sx)
    nnew <- nrow(xnew)
    my <- object$my
    sy <- object$sy
    
    lowrank <- object$lowrank
    hyper <- object$hyper
    
    oo <- .C("airGP_pred", xvar = as.double(x), yvar = as.double(obs), xnew = as.double(xnew), dim = as.integer(c(n, p, ncomp, lowrank$max.rank, nnew, nss)), hpar = as.double(hyper$sig), lpen = as.double(lowrank$lpenalty), ix.store = as.character(t(ix.store[ss,])), par.store = as.double(t(pars.store[ss,])), tolchol = as.double(lowrank$tol), fmeanstore = double((n + nnew) * nss), fstore = double((n + nnew) * nss), sigma.store = double(nss))
    
    return(list(f.hat = my + sy * matrix(oo$fmean, ncol = nss),
    f.samp = my + sy * matrix(oo$fstore, ncol = nss),
    sigma = sy * oo$sigma.store, subsamp = ss))
}

summary.airGP <- function(object, burn = 0.1, nmc = 200, ngraph = 10, use.names = FALSE, pars.plot = TRUE, graph.plot = TRUE, inclusion.plot = TRUE, keep.layout = FALSE, reorder = TRUE, ...){
    
    ix.store <- object$ix.store
    pars.store <- object$pars.store
    sigma.store <- object$sigma.store
    
    lamsq.grid <- object$lamsq.grid
    prox <- exp(-lamsq.grid * 0.01)
    nsamp <- nrow(pars.store)
    ncomp <- ncol(ix.store)
    
    n.postburn <- floor(nsamp * (1 - burn))
    nmc <- min(nmc, n.postburn)
    thin <- floor(n.postburn / nmc)
    ss <- nsamp - thin * ((nmc - 1):0)
    nss <- length(ss)
    
    nplots <- 2*pars.plot + graph.plot + inclusion.plot
    if(nplots == 0) keep.layout <- TRUE
    
    if(!keep.layout){
        plotpar <- par()
        nrow <- floor(sqrt(nplots))
        ncol <- ceiling(nplots/nrow)
        par(mfrow = c(nrow,ncol), font.main = 1)
    }
    if(pars.plot){
        
        lamsq.samp <- pars.store[ss,2*(1:ncomp),drop = FALSE]
        prox.samp <- exp(-lamsq.samp * 0.01)
        rhosq.samp <- pars.store[ss,2*(1:ncomp) - 1,drop = FALSE]
        Rsq.samp <- fnA(rhosq.samp)
        
        par.ix <- 1:ncomp
        if(reorder) par.ix <- order(apply(Rsq.samp, 2, quantile, pr = .75), decreasing = TRUE)
        boxplot(prox.samp[,par.ix], outline = FALSE, col = 3, ylim = c(0,1), ann = FALSE, bty = "n", names = par.ix, las = 3, cex.axis = 0.5)
        title(xlab = "Component", ylab = "Posterior spread", main = "prox")
        abline(h = prox, lty = 1, col = 2)
        
        boxplot(Rsq.samp[,par.ix], outline = FALSE, col = 3, ann = FALSE, bty = "n", ylim = c(0,1), names = par.ix, las = 3, cex.axis = 0.5)
        title(xlab = "Component", ylab = "Posterior spread", main = "Rsq")
        #abline(h = Rsq, lty = 1, col = 2)
    }
    
    p <- ncol(object$x)
    ngraph <- min(ngraph, p)
    gp.cocc <- matrix(0,p,p)
    v.incl <- rep(0,p)
    for(jj in ss){
        jj.incl <- rep(FALSE, p)
        for(k in 1:ncomp){
            if(ix.store[jj, k] != ""){
                ixn <- uncollapse(ix.store[jj, k], collapse = "\\.", mode = "numeric")
                if(max(ixn) > p) print(ix.store[jj,k])
                if(length(ixn) > 0){
                    gp.cocc[ixn,ixn] <- gp.cocc[ixn,ixn] + fnA(pars.store[jj, 2*(k-1)+1]) ##1
                    jj.incl[ixn] <- pmax(jj.incl[ixn], fnA(pars.store[jj, 2*(k-1)+1])) ##TRUE
                }
            }
        }
        v.incl <- v.incl + jj.incl
    }
    v.incl <- v.incl / length(ss)
    if(use.names){
        dimnames(gp.cocc)[[1]] <- dimnames(gp.cocc)[[2]] <- dimnames(object$x)[[2]]
    } else {
        dimnames(gp.cocc)[[1]] <- dimnames(gp.cocc)[[2]] <- 1:p
    }
    gp.cocc <- gp.cocc / length(ss)
    gp.cocc <- 20 *  gp.cocc / max(gp.cocc)
    
    if(inclusion.plot){
        plot(v.incl, ty = "h", bty = "n", xlab = "Predictor index", ylab = "Inclusion probability", main = "Predictor importance", ylim = c(0,1.1), xlim = c(0.5,p+0.5), col = tcol("darkgreen", .5))
        iz <- order(v.incl, decreasing = TRUE)[1:ngraph]
        if(use.names) text(iz, v.incl[iz] + 0.05, dimnames(object$x)[[2]][iz], pos = 4, cex = 1 - sqrt(nplots)*0.2, srt = 90, offset = 0)
    }
    
    iz <- sort(order(diag(gp.cocc), decreasing = TRUE)[1:ngraph])
    iz.conn <- unique(unlist(apply(gp.cocc[iz,] > 2, 1, which)))
    iz <- unique(c(iz, iz.conn))
    nz <- length(iz)
    gp.cocc.new <- gp.cocc[iz,iz,drop = FALSE]
    if(ngraph < nz) gp.cocc.new[ngraph + 1:(nz-ngraph), ngraph + 1:(nz-ngraph)] <- 0
    graph <- graph_from_adjacency_matrix(gp.cocc.new, mode = "undirected", weighted = TRUE, diag = FALSE)
    lo <- layout.fruchterman.reingold(graph, dim=2)
    if(graph.plot){
        plot(graph, layout=lo, vertex.size=diag(gp.cocc[iz,iz,drop = FALSE]), edge.width=E(graph)$weight, vertex.color = rep(c("orange","lightblue"), c(ngraph, nz - ngraph)), vertex.frame.color = NA, vertex.label.cex = (diag(gp.cocc[iz,iz,drop = FALSE])/20)^.3)
        title(main = "Predictor interaction")
    }
    
    if(!keep.layout) suppressWarnings(par(plotpar))
    invisible(list(graph = gp.cocc, inclusion = v.incl))
}

tcol <- Vectorize(function(col, alpha = 1) {x <- col2rgb(col)/255; return(rgb(x[1],x[2],x[3],alpha = alpha))}, "col")

## Utility function
uncollapse <- function(str, collapse = "", mode = "character"){
    a <- unlist(strsplit(str, collapse))
    mode(a) <- mode
    return(a)
}
fnA <- function(u) return(u/(1 + u))
fnB <- function(a, b) {
    ##vv <- pmax((1 - (1 - a)/(1 - b)), 0)
    vv <- exp(-9 * (1-a)/(1-b))
    return(vv)
}
lengthscale.fixation <- function(hyper){
    prox <- hyper$prox
    if(is.null(prox)) prox <- c(0.99, 0.93, 0.75, 0.45, 0.15)
    lamsq <- -log(prox) / 0.01
    pplam <- -diff(pbeta(c(1, prox, 0), hyper$lam[1], hyper$lam[2]), lag = 2)
    lplam <- log(pplam / sum(pplam))
    return(list(lamsq = lamsq, lplam = lplam))
}
magnitude.fixation <- function(hyper){
    rho.b <- hyper$rho[1]; rho.s <- rho.b*(1 - hyper$ERsq)/hyper$ERsq; rho.k <- ceiling(hyper$rho[2]); rho.del <- hyper$delta
    pmix <- function(rhosq) return(0.5 * pgamma(rhosq, rho.b, rho.s) + 0.5 * pgamma(rhosq, rho.b/rho.k, rho.s))
    my.eps <- 1e-6/rho.k
    #tau.grid <- unique(c(seq(rho.del, 1 - rho.del/2, rho.del), pmax(1 - rho.del, 1 - c(1e-1, 1e-2, 1e-4, 1e-6)/rho.k)))
    dseq <- seq(rho.del,1 - rho.del, len = 1 + ceiling((1/rho.del - 2)))
    #tau.max <- 1 + log(1 - rho.del)/rho.k
    #rho.M <- ceiling(log(-rho.k*rho.del / log(1 - rho.del), 2))
    rho.M <- ceiling(log(rho.del / 0.01, 2))
    tail <- rho.del/2^(1:rho.M)
    tau.grid <- unique(sort(c(tail, dseq, 1 - tail)))
    tau.grid <- tau.grid[tau.grid > pmix(10*my.eps)]
    tau.grid.mid <- 0.5 * (tau.grid[-1] + tau.grid[-length(tau.grid)])
    rhosq <- c(my.eps, sapply(tau.grid, function(u) uniroot(function(x) (pmix(x) - u), interval = c(my.eps/10, 1000), tol = 1e-10)$root))
    rhosq.cap <- c(0, 10*my.eps, sapply(tau.grid.mid, function(u) uniroot(function(x) (pmix(x) - u), interval = c(my.eps/10, 1000), tol = 1e-10)$root), Inf)
    pprho <- diff(c(pgamma(rhosq.cap, rho.b/rho.k, rho.s)))
    lprho <- log(pprho)
    lprho <- lprho - logsum(lprho)
    return(list(rhosq = rhosq, lprho = lprho))
}
logdiff <- function(lx) return(lx[-length(lx)] + log(expm1(diff(lx))))
logsum <- function(lx) return(max(lx) + log(sum(exp(lx - max(lx)))))
beta.pars <- function(target, p = c(0.5, 0.9)){
    f <- function(par) return(sum((target - qbeta(p, exp(par[1]), exp(par[2])))^2))
    f.opt <- optim(c(0,0), f, method = "BFGS")
    return(list(par = exp(f.opt$par), obj = f.opt$value))
}

