library(devtools)
#install_github("petrkeil/spasm")
library(spasm)
library(tidyverse)
library(vegan)
library(plyr)
library(foreach)
library(doMC)
library(ggrepel)



# ------------------------------------------------------------------------------
step_C_sim2 <- function(m)
{
  t(apply(X = m, MARGIN = 1, FUN = sample))
}


pk.dist <- function(m, form)
{
  x <- m[1,] # species x
  y <- m[2,] # species y
  
  
  # matching components
  a = sum(colSums(m) == 2)
  b = sum(x * (colSums(m) == 1) )
  c = sum(y * (colSums(m) == 1) )
  d = sum(colSums(m) == 0)
  n = a + b + c + d

  D = eval(parse(text = form))
  return(D)
}


# function that calculates the Z-score
Z_score <- function(m, algorithm, N.sim, form)
{
  res <- array(dim=c(nrow(m), nrow(m), N.sim))
  
  for(i in 1:N.sim)
  {
    # radnomize using a given algorithm
    null.m <- do.call(algorithm, list(m))
    # calculate the metric
    null.metric <- as.matrix(pk.dist(null.m, form))
    res[,,i] <- null.metric
  }
  
  # calculate the Z-score
  obs <-pk.dist(m, form)
  mean.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = mean))
  sd.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = sd))
  
  # the Z-score
  Z <- (obs - mean.null) / sd.null
  return(c(Z=Z, obs=obs))
}


# function that calculates the Z-score
Z_score_old <- function(m, algorithm, N.sim, form)
{
  res <- array(dim=c(nrow(m), nrow(m), N.sim))

  for(i in 1:N.sim)
  {
    # radnomize using a given algorithm
    null.m <- do.call(algorithm, list(m))
    # calculate the metric
    null.metric <- as.matrix(vegan::designdist(null.m,
                                               method = form,
                                               abcd = TRUE,
                                               terms = "binary"))
    res[,,i] <- null.metric
  }

  # calculate the Z-score
  obs <-vegan::designdist(m,
                          method = form,
                          abcd = TRUE,
                          terms = "binary")
  mean.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = mean))
  sd.null <- as.dist(apply(X = res, MARGIN = c(1,2), FUN = sd))

  # the Z-score
  Z <- (obs - mean.null) / sd.null
  return(c(Z=Z, obs=obs))
}

#m <- matrix(c(1,1,1,0,0,
#              0,0,1,1,0), byrow=TRUE, nrow=2, ncol=5)
#Z_score_old(m, algorithm="step_C_sim2", N.sim=1000, form = "a/n")
#Z_score(m, algorithm="step_C_sim2", N.sim=1000, form = "a/n")


################################################################################
# FORMULAS
################################################################################

list.of.formulas <- read.csv("../list_of_indices.csv") 
formulas <- select(list.of.formulas, formula, type, index_no)
formulas$formula <- as.character(formulas$formula)

# replace all 'n's with '(a+b+c+d)'
gsub(x = formulas$formula, pattern = "n", replacement = "(a+b+c+d)") 

################################################################################ 

params <- expand.grid(var.consp   = c(0.001, 0.01),
                      alpha = seq(-20, 20, by=5),
                      grain = c(16, 8),
                      N1 = c(100),
                      N2 = c(100))
# params <- split(params, 1:nrow(params))

registerDoMC(cores = 2)

output <- list()
for (i in 1:nrow(formulas))
{
  print(i)
  form <- as.character(formulas$formula[i])
  message(form)
  
  output[[i]] <- foreach(j = 1:nrow(params), .combine = "rbind") %dopar%
  {
      #message(paste("Simulation", j, "out of", nrow(params)))

      a <- sim.pair(abund.vect = c(params$N1[j], params$N2[j]),
                    var.consp = params$var.consp[j],
                    alpha = params$alpha[j],
                    plot.comm = FALSE)

      n.cells.side <- sqrt(params[j, 'grain'])

      m <- spasm::ppp.to.comm(a, dim.yx = c(n.cells.side, n.cells.side))$abundance

      # convert to binary matrix
      m.bin <-  m
      m.bin[m >= 1] <- 1

      Z <- Z_score2(m.bin, algorithm = "step_C_sim2", N.sim = 50, form = form)[]
      #Z <- Z_score_pair_analytic(m.bin, form = form)[]
      
      data.frame(Z = Z['Z'], obs = Z['obs'],
                             params[j,],
                             formula = form,
                             type = formulas$type[i])
    }

}


res <- ldply(output)
write.csv(res, file = "../results/simulation_results.csv", row.names=FALSE)
# ------------------------------------------------------------------------------

res <- read.csv("../results/simulation_results.csv")

res <- na.omit(res)


for.plot <- ddply(.data = res,
                  .variables = c("formula", "type"),
                  .fun=summarise,
                  Z_score = cor(alpha, Z, method= "kendall"),
                  Raw_metric = cor(alpha, obs, method = "kendall")) %>%
            gather(key = "raw_or_Z", value = "value", Z_score, Raw_metric )


match <- filter(.data = for.plot,
                type == "matching components")


png("../results/simulation_results.png", width=1200, height=700, res=190)
  ggplot(data=for.plot, aes(x = type, y = abs(value))) + 
   # geom_violin(aes(fill = color), scale = "width") +
    geom_boxplot(width = 0.5, aes(fill = type)) +
   # stat_summary(fun.y=median, geom="point", size=2, aes(colour = type)) +
    geom_label_repel(data = match, aes(x = type, y = abs(value), label = formula), 
               nudge_x = 0.1) +
    facet_grid(cols = vars(raw_or_Z)) +
    theme_bw() + #coord_flip() +
    theme(legend.position="none") +
    labs(y = "|Kendall's rank correlation with truth|", x = "")
dev.off()
