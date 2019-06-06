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


# function that calculates the Z-score
Z_score2 <- function(m, algorithm, N.sim, form)
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


################################################################################
# FORMULAS
################################################################################

formulas <- read.csv("../list_of_indices.csv") 


# REAL FORMULAS ----------------------------------------------------------------
# the actual indices:
real.formulas <- c(
# FROM KOLEFF ET AL (2003)
"(a+b+c)/((2*a+b+c)/2)" , # Whittaker (1962)
"(b+c)/2" , # Cody (1975), Lande (1996)
"(b+c)" , # Weiher & Boylen (1994)
"((a+b+c)^2)/((a+b+c)^2 - 2*b*c)", # Routledge (1977)
"(b+c)/(2*a + b +c)", # Wilson & Shmida (1984)
"(2*a + b +c)*(1-(a/(a+b+c)))", # Magurran (1988)
"1-(a*(2*a + b + c)/(2*(a+b)*(a+c)))", # Cody (1993)
"(b+c)/(a+b+c)", # Colwell & Coddington (1994)
"min(b,c)/(a+b+c)", # Williams (1996)
"(b*c + 1)/(((a+b+c)^2-(a+b+c))/2)", # Williams (1996)
"a/(a+c)", # Ruggiero et al. (1988)
"min(b,c)/(min(b,c) + a)", # Lennon et al. (2001)
"(2*abs(b-c))/(2*a + b + c)", # Lennon et al. (2001)
"1-( log((2*a + b + c)/(a + b + c))/log(2) )", # Harte & Kinzig (1997)
# FROM COOCCURRENCE ANALYSIS
"b*c", # classical C-score, Stone & Roberts (1990)
"a*d", # Togetherness, Stone & Roberts (1992)
# FROM HUBALEK (1982):  
"a/(max(a + b, a + c))", # Braun-Blanquet (1932)
"a/(min(a + b, a + c))", # Simpson (1943)
"a/(b + c)", # Kulczynski (1927)
"a/(a + b + c)", # Jaccard (1901)
"a/(a + 0.5*(b+c))", # Dice (1945), Sorensen (1948)
"a/(a + 2*(b + c))", # Sokal & Sneath (1963)
"0.5*(a/(a+b) + a/(a+c))", # Kulczynski (1927)
"(a/2)*(1/(a+b) + 1/(a+c))", # Driver & Kroeber (1932)
"a/(a+b) + a/(a+c)", # Johnson (1967)
"(a*a - b*c)/((a+b)*(a+c))", # McConnaughey (1964)
"a / (((a+b)*(a+c))^0.5 )", # Driver & Kroeber (1932), Ochiai (1957)
"(a*a)/((a+b)*(a+c))", # Sorgenfrei (1959)
"a/(((a+b)*(a+c))^0.5) - 0.5*(max(a+b, a+c))", # Fager & McGowan (1963)
"a/(a + b + c + d)", # Russell & Rao (1940)
"a/(0.5*(a*b + a*c) + b*c)", # Mountford (1962)
"(a*b + b*c) / (a*b + 2*b*c + c*d)", # Peirce (1884)
"(a - (a+b)*(a+c)) / ((a+b)*(c+d)*(a+c)*(b+d))" ,# Eyraud (1936)
"0.25*(a/(a+b) + a/(a+c) + d/(c+d) + d/(b+d))", # Sokal & Sneath (1963)
"(a+d) / (b+c)", # Sokal & Sneath (1963)
"(a + d) / (a+b+c+d)", # Sokal & Michener (1958)
"(50*pi)^(-1)* asin(((a + d) / (a+b+c+d))^0.5)", # Goodall (1967), Austin & Colwell (1977)
"(a+d)/(a + 0.5*(b+c) + d)", # Sokal & Sneath (1963)
"(a+d)/(a + 2*(b+c+d))", # Rogers & Tanimoto (1960)
"(a+d-b-c)/(a+b+c+d)", # Hamann (1961)
"a*d/(((a+b)*(c+d)*(a+c)*(b+d))^0.5)", # Sokal & Sneath (1963)
"(a*d - b*c)/((a+c)*(b+d))", # Pierce (1884)
"((a+b+c+d)*((a*d - b*c)^2)/((a+b)*(c+d)*(a+c)*(b+d)))", # Pearson (1905)
"((a+b+c+d)*((a*d - b*c)^2)/((a+b)*(c+d)*(a+c)*(b+d))) /(((a+b+c+d) +  ((a+b+c+d)*((a*d - b*c)^2)/((a+b)*(c+d)*(a+c)*(b+d))) )^0.5)", # Pearson (1905)
"sqrt(2) * (a*d - b*c) / (((a*d - b*c)^2 - (a+b)*(c+d)*(a+c)*(b+d))^0.5)", # cole (1949)
" (a*d + b*c) / (((a+b)*(c+d)*(a+c)*(b+d))^0.5)", # Yule (1912), Pearson & Heron (1913)
"(a*d - b*c)^2 / ((a+b)*(c+d)*(a+c)*(b+d))" ,# Doolittle (1885), Pearson (1926)
"(sqrt(a*d) + a)/(sqrt(a*d) + a + b + c)", # Baroni-Urbani & buser (1976)
"(sqrt(a*d) + a - b - c)/(sqrt(a*d) + a + b + c)",
"(a*d - b*c) / (a*d + b*c)", # Yule (1900)
"(sqrt(a*d) - sqrt(b*c)) / (sqrt(a*d) + sqrt(b*c))", # Yule (1900)
"cos(180*sqrt(b*c) / (sqrt(a*d) + sqrt(b*c)))", # Pearson & Heron (1913)
"4*(a*d - b*c)/( (a+d)^2 + (b+c)^2 )", # Michael (1920)
"(a+b+c+d)*a/((a+b)*(a+c))", # Forbes(1907)
"log(a) - log((a+b+c+d)) - log((a+b)/(a+b+c+d)) - log((a+c)/(a+b+c+d))", # Gilbert & Wells (1966)
"((a+b+c+d)*a -(a+b)*(a+c)) / ((a+b+c+d) * min(a+b, a+c) - (a+b)*(a+c))", # Forbes (1925)
"((a+b+c+d)*a -(a+b)*(a+c)) / ((a+b+c+d)*a + (a+b)*(a+c))") # Tarwid (1960))



real.formulas <- data.frame(formula = real.formulas, type="existing index")

# ------------------------------------------------------------------------------
# put all formulas together in one data frame
formulas <- rbind(one.letter, real.formulas)


################################################################################ 

params <- expand.grid(var.consp   = c(0.001, 0.01, 0.1),
                      alpha = seq(-20, 20, by=2.5),
                      grain = c(32, 16, 8),
                      N1 = c(100, 1000),
                      N2 = c(100, 1000))
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

      Z <- Z_score2(m.bin, algorithm = "step_C_sim2", N.sim = 100, form = form)[]
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
