library(Rcpp)
n = 800
p = 1000
source("R/ADMM_proj.R")
hM.time = Sys.time()
for (i in 1:5) {
  M <- matrix(rnorm(n*p), nrow = n)
  Sigma = cov(M)
  hM1 = ADMM_proj(Sigma)
}
hM.time = Sys.time() - hM.time
print(hM.time)
# M <- matrix(rnorm(1000), nrow = 10)
# Sigma = cov(M)

# # 1

# hM1 = ADMM_proj(M)$mat
# hM.time = ADMM_proj(M)$time
# hM.time
# # 2

# hM2 = ADMM_proj_Rcpp(M)$mat
# hM.time2 = ADMM_proj_Rcpp(M)$time
# hM.time2

# print(norm(hM1-hM2, type = "F"))
# print(norm(hM1, type = "F"))
# print(norm(hM2, type = "F"))