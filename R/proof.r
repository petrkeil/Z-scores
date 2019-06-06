
m.segr <- matrix(c(1,1,0,0,
                   0,0,1,1), byrow=TRUE, nrow=2, ncol=4)

m.agrg <- matrix(c(1,1,0,0,
                   1,1,0,0), byrow=TRUE, nrow=2, ncol=4)

m.part <- matrix(c(1,1,0,0,
                   0,1,1,0), byrow=TRUE, nrow=2, ncol=4)

# all possible combinations for each species
all <- combn(4,2)
# bind them to a character string
x.all <- apply(X = all, MARGIN = 2, FUN = paste, collapse = " ")
y.all <- apply(X = all, MARGIN = 2, FUN = paste, collapse = " ")
xy.all <- expand.grid(x = x.all, y=y.all, stringsAsFactors = FALSE)

x.all <- data.frame(do.call(rbind, strsplit(xy.all$x, split = " ")))
names(x.all) <- c("x1", "x2")
y.all <- data.frame(do.call(rbind, strsplit(xy.all$y, split = " ")))
names(y.all) <- c("y1", "y2")

xy.all <- data.frame(x.all, y.all)
