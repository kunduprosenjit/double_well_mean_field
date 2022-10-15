## Code by Neil MacLaren, 2/24/2022
## Corrections 10/13/2022. Old code commented with `####`. 

library(igraph)

dw <- function(x, r) {
    "This is the basic double well differential equation. For finding minima/maxima."
    -(x - r[1])*(x - r[2])*(x - r[3])
}

double_well_coupled <- function(x, r1, r2, r3, D, A, dt, noise = NULL, stress = rep(0, length(x))) {
    "Calculates the next step in a double well simulation and network-based coupling. Variables are generally as for `double_well`, with `x` a vector of current states, D as the coupling strength, and A as the adjacency matrix. The `noise` argument if not NULL should be a function that generates a vector of random values of length equal to the length of x."
    if(is.null(noise)) {# rowSums; test for re-run
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + stress)*dt
    } else {
        deltax <- (-(x - r1)*(x - r2)*(x - r3) + D*colSums(A*x) + stress)*dt + noise
    }
    
    nextx <- x + deltax
    nextx
}

noise <- function(n, s, f = rnorm) {
    "A function to generate `n` random noise values from a distribution with standard deviation `s`. Currently set up only for Gaussian noise, but could be expanded to support more distributions."
    f(n, 0, s)
}

lowerstate <- function(x, cutoff = 2.3) {
    "A counting function: return 1 if the node is in the lower state and zero otherwise. The current default cutoff is set to 2.3, which is /approximately/ the value of x_i at which the node state will be 'pulled' towards the separatrix when r = (1, 4, 7), the current default."
    ifelse(x <= cutoff, 1, 0)
}

save_plots <- FALSE # TRUE
use_noise <- FALSE # TRUE
run_simulation <- FALSE # TRUE
x_is_kdiff <- FALSE # TRUE

outfile <- "./data/kdiff-D-results.rda"

dataloc <- "./data/networks/"
ext <- ".rda"
choices <- gsub(ext, "", list.files(dataloc)) # list of empirical networks
datafiles <- paste0(dataloc, list.files(dataloc))
for(datafile in datafiles) load(datafile)

stepsize <- 1e-2 # for incrementing D

if(run_simulation) {
    results <- list()
    for(i in 1:length(choices)) {
        choice <- choices[i]
        print(choice)
        g <- get(choice)
        A <- as_adj(g, type = "both", sparse = FALSE) # this eliminates edge weights -- only returns adjacency
        
        nnodes <- vcount(g)
        nedges <- ecount(g)
        k <- degree(g)
        CV <- sd(k)/mean(k)

        k_diff <- quantile(k, probs = 0.95) - quantile(k, probs = 0.05)

        ## Node Systems
        r <- c(1, 2, 5) # double well parameters
        if(use_noise) s <- 0.005 else s <- 0 # sd of noise process
        D <- 0.01 # this is the initial value of D
        p <- if(s > 0) 3*s else 0.015 # perturbation strength
        u <- rep(0, nnodes) # stress vector

        τ <- 75 # duration in time units
        dt <- 0.01 # step size for integration
        T <- τ/dt

        initialx <- rep(.1, nnodes) # initial values of x_i

        ## Simulation Parameters
        x <- initialx # an updating vector of x_i

        cutoff <- .1*nnodes # to stop simulation; production value is .1
        state_cutoff <- 1.5
        in_lowerstate <- V(g) # initial vector of nodes in the lower state (all of them)

        ## Storage Vectors
        n_lowerstate <- numeric() # storage vector for the number of nodes in the lower state
        Ds <- numeric() # storage vector for the state of the bifurcation parameter

        ## Integrate and Store
        j <- 1
        while(length(in_lowerstate) > cutoff) {
            X <- matrix(0, nrow = T, ncol = nnodes) # storage matrix for the integration step states of x
            
            ## Integration
            for(t in 1:T) {
                X[t, ] <- x
                x <- double_well_coupled(x, r[1], r[2], r[3], D, A, dt, noise(nnodes, s), u)
            }
            
            ## Exit condition
            in_lowerstate <- V(g)[which(lowerstate(x, cutoff = state_cutoff) == 1)]
            if(length(in_lowerstate) <= cutoff) break

            ## Store
            n_lowerstate[j] <- length(in_lowerstate)
            Ds[j] <- D

            ## Iterate
            D <- D + stepsize
            j <- j + 1
        }

        ## Make Results Data Frame
        df <- data.frame(graph = choice, CV = CV, kdiff = k_diff,
                         D = Ds, n_lowerstate = n_lowerstate)
        results[[i]] <- df
    }
    results <- do.call(rbind, results)
    save(results, file = outfile)
} else {
    load(outfile)
    ##stepsize <- results$D[2] - results$D[1]
}
    

parts <- split(results, factor(results$graph))

find_range <- function(df, stepsize = 1e-3) {
    nnodes <- df$n_lowerstate[1] # the initial D is below that at which any of these networks starts to show transitions
    firststate <- floor(nnodes - (nnodes*.1))
    first <- which(df$n_lowerstate < firststate)[1]
    if(is.na(first)) return(0)
    ## otherwise...
    last <- nrow(df)
    firstD <- df$D[first]
    lastD <- df$D[last] + stepsize
    lastD - firstD
}

dat <- lapply(parts, function(x) {
    c(rangeD = find_range(x, stepsize),
      CV = unique(x$CV), kdiff = unique(x$kdiff),
      nstages = length(unique(x$n_lowerstate)))
})

df <- as.data.frame(do.call(rbind, dat))
df$N <- sapply(rownames(df), function(x) vcount(get(x)))

                                        # Retrieve local min and max for doublewell function
#### r <- c(1, 4, 7) 
r <- c(1, 2, 5)
#### ymin <- optimize(dw, c(1, 7), r = r)$minimum 
ymin <- optimize(dw, c(1, 5), r = r)$minimum
#### ymax <- optimize(dw, c(1, 7), r = r, maximum = TRUE)$maximum
ymax <- optimize(dw, c(1, 5), r = r, maximum = TRUE)$maximum

df$k_min <- sapply(rownames(df), function(x) min(degree(get(x))))
df$k_max <- sapply(rownames(df), function(x) max(degree(get(x))))

df$xval <- ymin/df$k_min - ymax/df$k_max

                                        # alternate, y2/kmin - y1/kmax
ht <- 6; wd <- 6.5
if(save_plots) {
    pdf("./img/compare-kminmax-with-D.pdf", height = ht, width = wd)
} else dev.new(height = ht, width = wd)
par(mar = c(5, 6, 2, 2) + 0.1)
plot(rangeD ~ xval, data = df, pch = 1, col = "black", cex = 2, lwd = 1.75,
     axes = FALSE,
     #### xlim = c(0, 2.25),
     xlim = c(0, 1.5),
     ylim = range(df$rangeD),
     xlab = "", ylab = ""
     )
box()
axis(1, cex.axis = 1.75, las = 1)
axis(2, cex.axis = 1.75, las = 1)
title(xlab = expression(italic(tilde(y))^(2)/italic(k)[min] - italic(tilde(y))^(1)/italic(k)[max]),
      cex.lab = 1.75, line = 3.5)
title(ylab = "Magnitude of range of D", line = 4.5, cex.lab = 1.75)
mtext("(b)", line = -1.3, adj = .01, font = 1, cex = 1.5)
if(save_plots) dev.off()
