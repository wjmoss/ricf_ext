glist <- lapply(rep(20,100), generateModel, n=200, k=2, d=d, b=b)

cumt_o = 0
res_o = list()
for (i in 1:100){
  t = proc.time()[3]
  try(res_o[[i]] <- ricfR(L=t(glist[[i]]$B), O=glist[[i]]$Omega, X=t(glist[[i]]$Y)))
  cumt_o = cumt_o + proc.time()[3] - t
}


cumt_c = 0
res_c = list()
for (i in 1:100){
  t = proc.time()[3]
  try(res_c[[i]] <- ricf_R(L=t(glist[[i]]$B), O=glist[[i]]$Omega, S=glist[[i]]$Y %*% t(glist[[i]]$Y)))
  cumt_c = cumt_c + proc.time()[3] - t
}

cumt_s = 0
res_s = list()
for (i in 1:100){
  t = proc.time()[3]
  try(res_s[[i]] <- ricf_R_scc(L=t(glist[[i]]$B), O=glist[[i]]$Omega, S=glist[[i]]$Y %*% t(glist[[i]]$Y)))
  cumt_s = cumt_s + proc.time()[3] - t
}


cumt_a = 0
res_a = list()
for (i in 1:100){
  t = proc.time()[3]
  try(res_a[[i]] <- ricf_R_acc(L=t(glist[[i]]$B), O=glist[[i]]$Omega, S=glist[[i]]$Y %*% t(glist[[i]]$Y)))
  cumt_a = cumt_a + proc.time()[3] - t
}

sum(sapply(1:100, function(i) res_o[[i]]$converged))
sum(sapply(1:100, function(i) res_c[[i]]$converged))
sum(sapply(1:100, function(i) res_s[[i]]$converged))
sum(sapply(1:100, function(i) res_a[[i]]$converged))
