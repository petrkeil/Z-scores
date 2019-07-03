
# Code for the Z-score analysis
# Author: Petr Keil
# pkeil@seznam.cz


################################################################################ 
# R packages
################################################################################

# library(devtools)
# install_github("petrkeil/spasm")
library(spasm)
library(tidyverse)
library(vegan)
library(plyr)
library(foreach)
library(doMC)
library(ggrepel)
library(gridExtra)

################################################################################
# Randomization algorithm (null model).
################################################################################

step_C_sim2 <- function(m)
{
  t(apply(X = m, MARGIN = 1, FUN = sample))
}

################################################################################
# function that calculates the given index (so far works for 2 species only)
################################################################################

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

################################################################################
# function that calculates the Z-score
################################################################################

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


#m <- matrix(c(1,1,1,0,0,
#              0,0,1,1,0), byrow=TRUE, nrow=2, ncol=5)
#Z_score(m, algorithm="step_C_sim2", N.sim=1000, form = "a/n")


################################################################################
# Index formulas
################################################################################

list.of.formulas <- read.csv("../list_of_indices.csv") 
formulas <- select(list.of.formulas, formula, type, index_no)
formulas$formula <- as.character(formulas$formula)

# replace all 'n's with '(a+b+c+d)'
gsub(x = formulas$formula, pattern = "n", replacement = "(a+b+c+d)") 

################################################################################ 
# Simulation parameters
################################################################################

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      grain = c(32, 16, 8),
                      N1 = c(100, 1000),
                      N2 = c(100, 1000))

################################################################################ 
# Simulations
################################################################################

registerDoMC(cores = 4)

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

      Z <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = form)[]
      #Z <- Z_score_pair_analytic(m.bin, form = form)[]
      
      data.frame(Z = Z['Z'], obs = Z['obs'],
                             params[j,],
                             formula = form,
                             index_no = formulas$index_no[i],
                             type = formulas$type[i])
    }

}


res <- ldply(output)
write.csv(res, file = "../results/simulation_results.csv", row.names=FALSE)

################################################################################ 
# Plotting results
################################################################################

res <- read.csv("../results/simulation_results.csv")
res <- na.omit(res)

for.plot <- ddply(.data = res,
                  .variables = c("formula", "type"),
                  .fun=summarise,
                  Z_score = cor(alpha, Z, method= "spearman"),
                  Raw_metric = cor(alpha, obs, method = "spearman")) %>%
            gather(key = "raw_or_Z", value = "value", Z_score, Raw_metric )


match <- filter(.data = for.plot,
                type == "matching component ")



bx <- ggplot(data=for.plot, aes(x = type, y = abs(value))) + 
  # geom_violin(aes(fill = color), scale = "width") +
  geom_boxplot(width = 0.5, aes(fill = type), outlier.shape = 1) +
  # stat_summary(fun.y=median, geom="point", size=2, aes(colour = type)) +
  geom_label_repel(data = match, aes(x = type, y = abs(value), label = formula), 
                   nudge_x = 0.1) +
  facet_grid(cols = vars(raw_or_Z)) +
  theme_bw() + #coord_flip() +
  theme(legend.position="none") +
  labs(y = "| Spearman correlation with alpha |", x = "")

png("../results/simulation_results.png", width=1200, height=700, res=190)
  print(bx)
dev.off()

# boxplots medians
ddply(for.plot, .(type, raw_or_Z), summarise, med = median(abs(value)))

# boxplots inter-quartile range
ddply(for.plot, .(type, raw_or_Z), summarise, med = IQR(abs(value)))


# which are the extremely poorly performing indices?
x <- for.plot[for.plot$type == "existing index",]
x <- x[order(abs(x$value)),]
head(x)
x <- x[order(abs(x$value), decreasing = TRUE),]
head(x)


# which are the extremely well performing raw indices?
x <- for.plot[for.plot$raw_or_Z == "Raw_metric",]
x <- x[order(abs(x$value)),]
head(x)
x <- x[order(abs(x$value), decreasing = FALSE),]
head(x)



     



################################################################################ 
# Simulation parameters
################################################################################

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1, 1),
                      alpha = seq(-20, 20, by=2.5),
                      grain = c(32, 16, 8),
                      N1 = c(100, 1000),
                      N2 = c(100, 1000))

################################################################################ 
# Simulations
################################################################################

registerDoMC(cores = 4)


JACC <- foreach(j = 1:nrow(params), .combine = "rbind") %dopar%
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
    
    Z.a <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = "a")[]
    Z.b <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = "b")[]
    Z.c <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = "c")[]
    Z.d <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = "d")[]
    Z.jacc <- Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 400, form = "a / (a + b + c)")[]
    #Z <- Z_score_pair_analytic(m.bin, form = form)[]
    
    data.frame(a = Z.a['obs'],
               b = Z.b['obs'],
               c = Z.c['obs'],
               d = Z.d['obs'],
               Zscore_a = Z.a['Z'], 
               Zscore_b = Z.b['Z'],
               Zscore_c = Z.c['Z'], 
               Zscore_d = Z.d['Z'],
               Zscore_Jaccard = Z.jacc['Z'], 
               Jaccard = Z.jacc['obs'],
               params[j,])
}


plot(JACC$alpha, JACC$a)
cor(JACC$alpha, JACC$b, method = "spearman")


a_vs_jacc <- ggplot(data = JACC, aes(x = Jaccard, Zscore_a)) +
             geom_point(shape = 1) + theme_bw() + ggtitle("a")
a_vs_jacc


alpha_vs_jacc <- ggplot(data = JACC, aes(x = alpha, Jaccard)) +
  geom_point(shape = 1) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("b")
alpha_vs_jacc


alpha_vs_Z <- ggplot(data = JACC, aes(x = alpha, Zscore_a)) +
  geom_point(shape = 1) + geom_smooth(method = "lm", se = FALSE) + theme_bw() + ggtitle("c")
alpha_vs_Z


png("../results/Z_vs_Jaccard.png", width= 1200, height = 400, res=100)
grid.arrange(a_vs_jacc, alpha_vs_jacc, alpha_vs_Z, ncol=3)
dev.off()





a <- sim.pair(abund.vect = c(100, 100),
                        var.consp = 0.001,
                        alpha = 20,
                        plot.comm = TRUE)
m <- spasm::ppp.to.comm(a, dim.yx = c(8,8))$abundance
m.bin <-  m
m.bin[m >= 1] <- 1

# b and c components
sum(colSums(m.bin) == 2) # a
sum(colSums(m.bin) == 1 & m.bin[1,] == 1) # b
sum(colSums(m.bin) == 1 & m.bin[2,] == 1) # c
sum(colSums(m.bin) == 0) # d

m.bin
Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 200, form = "a")[]
Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 200, form = "b")[]
Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 200, form = "c")[]
Z_score(m.bin, algorithm = "step_C_sim2", N.sim = 200, form = "d")[]
