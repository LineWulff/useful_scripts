## input list of 5 sets of genes!! (you can repeat ine list several times if you don't have 5)

vennover <- function(x){
A <- x[[1]]
B <- x[[2]]
C <- x[[3]]
D <- x[[4]]
E <- x[[5]]
n12 <- intersect(A, B)
n13 <- intersect(A, C)
n14 <- intersect(A, D)
n15 <- intersect(A, E)
n23 <- intersect(B, C)
n24 <- intersect(B, D)
n25 <- intersect(B, E)
n34 <- intersect(C, D)
n35 <- intersect(C, E)
n45 <- intersect(D, E)
n123 <- intersect(n12, C)
n124 <- intersect(n12, D)
n125 <- intersect(n12, E)
n134 <- intersect(n13, D)
n135 <- intersect(n13, E)
n145 <- intersect(n14, E)
n234 <- intersect(n23, D)
n235 <- intersect(n23, E)
n245 <- intersect(n24, E)
n345 <- intersect(n34, E)
n1234 <- intersect(n123, D)
n1235 <- intersect(n123, E)
n1245 <- intersect(n124, E)
n1345 <- intersect(n134, E)
n2345 <- intersect(n234, E)
n12345 <- intersect(n1234, E)
return(n12345)
}