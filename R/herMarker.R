
#' Marker Based Heritability Calutation
#'
#' \code{herMarker} returns the marker based heritability using sommer package and RKHS function from BGLR package.
#'
#' @param geno The name of the file which the genotypic-data are to be read from or a matrix in the r-environment. Data coded as (-1,0,1)
#' @param samp The name of the file which the genotypes are to be read from or a vector  with the genotype names
#' @param phen The name of the file which the phenotypic-data are to be read from or a data.frame in the r-environment.
#' @param method A String c("RKHS","sommer")
#' @param traits NULL by default or a vector with trait names
#'
#' @return two components
#' @export
#'
#' @examples
#'
#' # library(tidyverse)
#' # library(sommer)
#' # data(DT_cpdata)
#'
#' # geno <- GT_cpdata
#' # samp <- rownames(GT_cpdata)
#' # phen <- DT_cpdata
#'
#' # tmp <- herMarker(geno, samp, phen, method = c("RKHS","sommer"), traits=NULL)
#'
#' # names(tmp$data)
#' #  tmp$data[,-c(4,5)] %>% spread(method,h)
#' #
#' # tmp$data %>%
#' #   ggplot(aes(x=method, y=h,fill=method, label=round(h,2)))+
#' #   geom_bar(stat = "identity", position = "dodge" )+
#' #   theme_bw()+
#' #   theme(axis.text.x = element_text(hjust = 1,angle = 75))+
#' #   geom_text(aes(method), size=3,nudge_y = 0.1)+
#' #   facet_wrap(~trait,ncol = 4,scales = "free_x")
#' @author Johan Aparicio, \email{j.aparicio@cgiar.org}

"herMarker" <- function(geno, samp, phen, method =  c("RKHS", "sommer"), traits=NULL){

  if (length(samp)>1) {
    samp <- samp
  } else {
    samp = read.delim(samp, header = F)[,1]
  }

  if (inherits(geno, what = "matrix")) {
    geno <- as.matrix(geno)
    G <- geno
    rownames(G) <- samp
  } else {
    geno <- geno
    G = read.table(geno, row.names = as.character(samp), header = F)
    G <- G[ order(row.names(G)) , ]
  }

  if (inherits(phen, what = "data.frame")) {
    phen <- as.data.frame(phen)
    genoname <<- names(phen)[1]
  } else {
    phen <- data.frame(read.csv(phen))
    genoname <<- names(phen)[1]
    phen <- arrange(phen, get(genoname))
  }

  mt <- c("RKHS", "sommer") # "ASReml"
  if (sum(method%in%mt)!=length(method)) {
    stop("Check the method names")
  }

  phen <- GS:::checkPhen(phen)

  # traits
  if(is.null(traits)) traits <- names(phen)[names(phen)!=genoname]

  traits <- intersect(traits, names(phen)[names(phen)!=genoname])

  message("[]==============================================================[]")
  message("[]======== Marker based heritability calculation ===============[]")
  message("[]=================  BGLR - sommer package =====================[]")
  message("[]======= Last update: 2020-03-22  Johan Aparicio ==============[]")
  message("[]==============================================================[]\n")

  Gen <- list()
  for (i in traits) {
    LinesA <- as.character(phen[!is.na(phen[,i]),genoname])  # A
    LinesB <- as.character(rownames(G))                      # B
    Gen[[i]] <- intersect(LinesA,LinesB)
    message(i, "\t == " ,length(intersect(LinesA,LinesB)))
  }


  # Dataframe with results --------------------------------------------------
  out_table = data.frame(trait = as.character(), method= as.character(),
                         h=as.numeric(), corr=as.numeric() ,  finishedAt = as.character())

  message("\ntrait","\tmethod", "\th","\tcorr", "\tfinishedAt" )
  #--------------------------------------------------------------------------

  modMar <- c("RKHS","sommer")
  stI_list <- matrix(data=list(), nrow=length(traits), ncol=length(modMar),
                     dimnames=list(traits, modMar))

  #--------------------------------------------------------------------------


  for (i in traits  ) {

    LinesA <- as.character(phen[!is.na(phen[,i]),genoname])  # A
    LinesB <- as.character(rownames(G))                      # B
    GT <- G[rownames(G)%in%intersect(LinesA,LinesB),]   # Marker

    phen2 <- droplevels(subset(phen,phen[,1]%in%rownames(GT)))
    phen2$var <- phen2[,i]
    phen2$level <- as.factor(phen2[,genoname])

    #------ ainverse  ---------


    G_hat <- GT
    AHAT<-sommer::A.mat(G_hat)
    AHAT_blend<-0.98*AHAT+0.02*diag(x=1,nrow=nrow(G_hat),ncol=nrow(G_hat))
    det(AHAT_blend)
    Ainv.blend<-solve(AHAT_blend)

    rownames(Ainv.blend) <- colnames(Ainv.blend) <- NULL
    AHAT.inv.sparse<-GS:::full2sparse(Ainv.blend)  #lower diag sparse matrix (Ainv.bend)
    colnames(AHAT.inv.sparse)<-c('Row','Column','Ainverse')

    # Preparing Ginverse
    ahatinv<-data.frame(AHAT.inv.sparse)
    attr(ahatinv,"rowNames")<-as.character(phen2[,genoname]) # A vector with ID names in order
    attr(ahatinv,"colNames")<-as.character(phen2[,genoname])
    attr(ahatinv,"INVERSE")<-TRUE  # Very important!

    phen2[,genoname] <- as.factor(phen2[,genoname])

    #-----------------------

    y <- phen2$var
    n <- length(y)

    h_rkhs <- h_som_1 <- h_som_2 <- c()
    time2  <- time4 <- time5 <- c()
    corrk2 <- corrk4 <- corrk5 <- c()
    fm2 <- fm3 <- fm4 <- list()


    # Models ------------------------------------------------------------------

    yNA<-y

    # ------ RKHS

    if("RKHS"%in%method){


      D = as.matrix( dist( GT, method="euclidean")) ^2
      D = D / mean(D)

      h = 0.5                                # fixed = list(~coord,data=phen2,model='FIXED')
      K = exp(-h * D)                        # fixed = list(~factor(Haplotype),data=phen2,model='FIXED')
      ETA = list(list(model="RKHS",   K=K))  # fixed = list(~factor(Meso.or.Andean),data=phen2,model='FIXED')


      fm2 =BGLR::BGLR( y=yNA, ETA=ETA,
                       nIter=10000, burnIn=1000, thin=5,
                       verbose=FALSE)

      h_rkhs <- fm2$ETA[[1]]$varU/(fm2$ETA[[1]]$varU + fm2$varE)
      h_rkhs <- round(h_rkhs,3)

      corrk2 <- cor(y,fm2$yHat) %>% round(.,3)
      time2 <- paste0(Sys.time())

      message(i, "\tRKHS","\t",h_rkhs,"\t",corrk2, "\t", time2 )

    }

    # ------ sommer 2

    if("sommer"%in%method){

      phen2$vSom <- yNA
      phen2$level2 <-phen2$level
      DHAT_blend <-sommer::D.mat(G_hat)

      fm4 = sommer::mmer(vSom~1,
                         random=~sommer::vs(level,Gu=AHAT_blend)+sommer::vs(level2,Gu=DHAT_blend),
                         rcov=~units,
                         data=phen2,verbose=FALSE )

      h_som_1 <- round(as.numeric(sommer::pin(fm4, h2 ~ (V1) / ( V1+V3) )[1]) ,3 )
      h_som_2 <- round(as.numeric(sommer::pin(fm4, h2 ~ (V1+V2) / ( V1+V2+V3) )[1]) ,3 )

      fm4$U$`u:level`$vSom <- as.data.frame(fm4$U$`u:level`$vSom)
      rownames(fm4$U$`u:level`$vSom) <- gsub("level","",rownames(fm4$U$`u:level`$vSom))

      corrk4 <- cor(fm4$U$`u:level`$vSom[,],y, use="complete") %>% round(.,3)
      corrk5 <- cor(fm4$U$`u:level2`$vSom,y, use="complete") %>% round(.,3)

      time4 <- paste0(Sys.time())
      time5 <- time4

      message(i,"\ts_narrow","\t" ,h_som_1, "\t",corrk4,"\t", time4)
      message(i,"\ts_broad", "\t" ,h_som_2, "\t",corrk5,"\t", time5)
    }


    # final results ---------

    hf <- c(h_rkhs,h_som_1,h_som_2)
    corf <- c(corrk2,corrk4,corrk5)
    time <- c(time2,time4,time5)

    if("sommer"%in%method) meth <- c(sort(method)[1],"s_narrow","s_broad")
    else meth  <- method

    tmp = data.frame(trait = i,
                     method=meth ,
                     h = hf,
                     corr = corf,
                     finishedAt = time)
    out_table <- rbind(out_table,tmp)

    stI_list[[i, 1]] <- fm2
    stI_list[[i, 2]] <- fm4

  }


  message("\n[]============================ End =============================[]\n")

  return(list(data = out_table,
              models = stI_list ))

}


