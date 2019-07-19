# setup({
#   nonsensepars <- list(foeH=0, foeL=0, pH = 0)
# })

test_that('rfoi returns correct object', {
  testlen <- 1057; tfoeH = .3; tfoeL = .1; tpH = .2
  rn <- rfoi(testlen, tfoeH, tfoeL, tpH)
  expect_equal(length(rn), testlen)
  expect_equal(all((rn == tfoeH) | (rn == tfoeL)), TRUE)
  whchH <- (rn == tfoeH)
  expect_equal(sum(whchH) < sum(!whchH), TRUE) # while possible this could fail, extremely unlikely
})

test_that('seros returns correctly sized, ordered object', {
  nCircs <- 10; nPop <- 10; mAge <- 10
  sr <- seros(nCircs, mAge, nPop)
  expect_equal(dim(sr)[2], mAge)
  expect_equal(dim(sr)[1], nCircs*nPop)
  expect_equal(all(sr[1,] == sr[1+nCircs,]), TRUE)
})

test_that('expfoi creates matching size, order object', {
  nCircs <- 10; nPop <- 10; mAge <- 10
  rn <- rfoi(nPop, .2, .1, .1)
  res <- expfoi(rn, nCircs, mAge)
  expect_equal(dim(res)[2], mAge)
  expect_equal(dim(res)[1], nCircs*nPop)
  expect_equal(length(unique(c(res[1:nCircs,]))), 1)
  expect_equal(unique(c(res[1:nCircs,])), rn[1])
})

test_that('pruned series have correct properties', {
  fake <- matrix(1:4, ncol=3, nrow=4)
  p1 <- prune.exposures(fake)
  p2 <- prune.exposures(cbind(0,fake))
  p3 <- prune.exposures(matrix(1:2, ncol=3, nrow=4))
  p4 <- prune.exposures(matrix(1:4, ncol=4, nrow=4, byrow = T))
  p5 <- prune.exposures(matrix(1:4, ncol=8, nrow=4, byrow = T))
  expect_equal(all(p1[,1]==1:4), TRUE)
  expect_equal(all(p1[,-1]==0), TRUE)
  expect_equal(all(p2[,2]==1:4), TRUE)
  expect_equal(all(p2[,-2]==0), TRUE)

  expect_equal(all(apply(unique(p4, MARGIN = 2), 1, length)==3), TRUE)
  expect_equal(all(p4!=2 & p4!=4), TRUE)
  expect_equal(all(apply(unique(p5, MARGIN = 2), 1, length)==5), TRUE)
})

# test_that('nPxA produces correctly dimensioned results', {
#   nCirc <- 1; nPop <- 1; mA <- 10
#   dt <- nPxA(nonsensepars, nCirculation = nCirc, nPop = nPop, maxAge = mA)
# })
