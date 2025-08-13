distance_MuFiMeshGP <- function(x1, x2 = NULL) {
  if (is.null(x2)) D <- anisoDist1_cpp(x1) else D <- anisoDist2_cpp(x1, x2)
  return(D)
}
