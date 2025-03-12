## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
# Load the package and dependencies
library(ppgm)

## ----load occurrence data-----------------------------------------------------
# load Sceloporus occurrence data
data(occurrences)

## ----view occurrence----------------------------------------------------------
# View the occurrence data to check all in order
head(occurrences,n=3)

## ----load trees---------------------------------------------------------------
# Load the Sceloporus phylogenies and select only 10
data(sampletrees)
vigtrees <- sample(sampletrees, size=10)

## ----load fossils-------------------------------------------------------------
data("scel_fossils")
head(scel_fossils)

## ----examine bio1-------------------------------------------------------------
bio1occ <- getBioclimVars(occurrences = occurrences, which.biovars=1)
summary(bio1occ$bio1)

## ----examine bio4-------------------------------------------------------------
bio4occ <- getBioclimVars(occurrences = occurrences, which.biovars=4)
summary(bio4occ$bio4)


## ----set BM bounds------------------------------------------------------------
# Set bounds for Brownian Motion
bmbounds <- list(sigsq = c(min = 0, max = 10000))

## ----set OU bounds------------------------------------------------------------
# Set bounds for OU
oubounds <- list(sigsq = c(min = 0, max = 1000000),alpha = c(min = 0, max = 100000))

## ----set control--------------------------------------------------------------
#set control for model
contr <- list(niter=20)

## ----run ppgm 1 BM no fossils-------------------------------------------------
ppgm_bm <- ppgm(occurrences, trees = vigtrees, bounds = bmbounds, control = contr, 
                    which.biovars = c(1,4), model = "BM",  path = "BM")


## ----image-ppgm-nofossils-bm-bio1, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram showing Mean Annual Temperature', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/BMbio1.jpg")

## ----image-ppgm-nofossils-bm-bio4, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram showing Temperature Seasonality', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/BMbio4.jpg")

## ----run ppgm 2 OU no fossils-------------------------------------------------
ppgm_ou <- ppgm(occurrences, trees = vigtrees, bounds = oubounds, control = contr, 
                    which.biovars = c(1,4), model = "OU", path = "OU")


## ----image-ppgm-nofossils-OU-bio1, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram showing Mean Annual Temperature', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/OUbio1.jpg")

## ----image-ppgm-nofossils-OU-bio4, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram showing Temperature Seasonality', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/OUbio4.jpg")

## ----climate envelope output--------------------------------------------------
head(ppgm_bm$cem)

## ----climate envelope output2-------------------------------------------------
head(ppgm_bm$envelope[[1]])

## -----------------------------------------------------------------------------
head(ppgm_bm$treedata_min[[1]]$data)

## ----format output for ppgmMESS-----------------------------------------------
cem_min <- cbind(ppgm_bm$cem[,1], ppgm_bm$cem[,2])
cem_max <- cbind(ppgm_bm$cem[,5], ppgm_bm$cem[,6])
rownames(cem_min) <- rownames(cem_max) <- rownames(ppgm_bm$cem)

## ----run ppgmMESS-------------------------------------------------------------
mess <- ppgmMESS(cem_min, cem_max, ppgm_bm$node_est, tree = vigtrees, timeslice = 15,
                 which.biovars = c(1,4), which.plot = "none")

## ----image-ppgm-15mess, echo = FALSE, message=FALSE, fig.align='center', fig.cap='MESS at 15mya', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/MESS15Multi.jpg")

## -----------------------------------------------------------------------------
vigfossils <- head(scel_fossils)

## ----ppgm run with fossils----------------------------------------------------
ppgm_bmfossil <- ppgm(occurrences, trees = vigtrees, fossils = vigfossils, bounds = bmbounds,
                      control = contr, which.biovars = c(1,4), model = "BM",
                      path = "BMfossil", permut = 5)

## ----image-ppgm-nofossils-bm-bio1-again, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram without fossils showing Mean Annual Temperature', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/BMbio1.jpg")

## ----image-ppgm-fossils-bm-bio1, echo = FALSE, message=FALSE, fig.align='center', fig.cap='Traitgram with fossils showing Mean Annual Temperature', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/BMfossilbio1.jpg")

## ----format fossil output for ppgmMESS----------------------------------------
fcem_min <- cbind(ppgm_bmfossil$cem[,1], ppgm_bmfossil$cem[,2])
fcem_max <- cbind(ppgm_bmfossil$cem[,5], ppgm_bmfossil$cem[,6])
rownames(fcem_min) <- rownames(fcem_max) <- rownames(ppgm_bmfossil$cem)

## ----run fossil ppgmMESS------------------------------------------------------
mess_fossil <- ppgmMESS(fcem_min, fcem_max, ppgm_bmfossil$node_est, tree=vigtrees, timeslice=15,
                        which.biovars = c(1,4), which.plot = "none", fossils = vigfossils,
                        path = "fossil")

## ----image-ppgm-15mess-nofossils, echo = FALSE, message=FALSE, fig.align='center', fig.cap='MESS at 15mya, no fossils', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/MESS15Multi.jpg")

## ----image-ppgm-15mess-fossils, echo = FALSE, message=FALSE, fig.align='center', fig.cap='MESS at 15mya, fossils', out.width="75%", fig.pos='H'----
knitr::include_graphics("images/fossilMESS15Multi.jpg")

