Ginverse <- function(G){
  # additive relationship matrix
  A <-sommer::A.mat(G)
  AHAT_blend<-0.98*A+0.02*diag(x=1,nrow=nrow(A),ncol=nrow(A))
  Ainv.blend<-solve(AHAT_blend)
  rownames(Ainv.blend) <- colnames(Ainv.blend) <- NULL
  AHAT.inv.sparse<-GS::full2sparse(Ainv.blend)
  colnames(AHAT.inv.sparse)<-c('Row','Column','Ainverse')

  # Preparing Ginverse
  ahatinv<-data.frame(AHAT.inv.sparse)
  attr(ahatinv,"rowNames")<-as.character(rownames(A)) # A vector with ID names in order
  attr(ahatinv,"colNames")<-as.character(colnames(A))
  attr(ahatinv,"INVERSE")<-TRUE
  return(ahatinv)
}
