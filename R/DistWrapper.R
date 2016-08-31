ddist_Uni <- function(y, Theta, Dist, log = FALSE){
  ddist_univ(y, Theta, Dist, log)
}

rdist_Uni <- function(Theta, Dist){
  rdist_univ(Theta, Dist)
}

pdist_Uni <- function(q, Theta, Dist){
  pdist_univ(q, Theta, Dist)
}

qdist_Uni <- function(p, Theta, Dist){
  pdist_univ(p, Theta, Dist)
}

mdist_Uni <- function(Theta, Dist){
  mdist_univ(Theta, Dist)
}

Score_Uni <- function(y, Theta, Dist){
  Score_univ(y, Theta, Dist)
}

IM_Uni <- function(Theta, Dist){
  IM_univ(Theta, Dist)
}


ddist_Multi <- function(y, Theta, Dist, log = FALSE){
  ddist_multi(y, Theta, length(y), Dist, log)
}

rdist_Multi <- function(Theta, N, Dist){
  rdist_multi(Theta, N, Dist)
}

Score_Multi <- function(y, Theta, Dist){
  Score_multi(y, Theta, length(y), Dist)
}
