manski2 <- function(data1,control.M1,control.M2,control.M2a,treated.M1,treated.M2,treated.M2a)  {

  temp <- rep(NA,11)

  # save true ATE estimate in this sample
  temp[5] <- mean(data1$Y1[data1$Tr == 1]) - mean(data1$Y0[data1$Tr == 0])


  # calculate lower and upper Manski bounds for Y0 (CONTROLS)
  denom.Co <- sum(data1$Tr == 0)

  w1.Co <- sum(data1$Tr == 0 & data1$M1 == 0) / denom.Co                    # weight for fully observed outcomes
  w2.Co <- sum(data1$Tr == 0 & data1$M1 == 1 & data1$M2 == 0) / denom.Co    # weight for outcomes observed during follow-up
  w3.Co <- sum(data1$Tr == 0 & data1$M1 == 1 & data1$M2 == 1) / denom.Co    # weight for outcomes imputed with extreme values

  # consistency check
  stopifnot(abs(w1.Co+w2.Co+w3.Co-1) < .00001)

  m1.Co <- mean(data1$Y[data1$Tr == 0 & data1$M1 == 0])
  m2.Co <- mean(data1$Y[control.M2a & data1$M2 == 0])

  # consistency check
  if(length(data1$Y[control.M2a & data1$M2 == 0]) == 0) m2.Co <- 0

  # lower bound
  temp[3] <- w1.Co*m1.Co + w2.Co*m2.Co + w3.Co*0

  # upper bound
  temp[4] <- w1.Co*m1.Co + w2.Co*m2.Co + w3.Co*1

  # calculate less naive estimator that uses randomization and follow-up outcomes but ignores 2nd stage missingness
  temp[7] <- w1.Co/(w1.Co+w2.Co)*m1.Co + w2.Co/(w1.Co+w2.Co)*m2.Co

  # calculate Manski bounds for naive estimator
  temp[8] <- m1.Co*w1.Co
  temp[9] <- m1.Co*w1.Co + (w2.Co+w3.Co)


  # calculate lower and upper Manski bounds for Y1 (TREATED)
  denom.Tr <- sum(data1$Tr == 1)

  w1.Tr <- sum(data1$Tr == 1 & data1$M1 == 0) / denom.Tr                    # weight for fully observed outcomes
  w2.Tr <- sum(data1$Tr == 1 & data1$M1 == 1 & data1$M2 == 0) / denom.Tr    # weight for outcomes observed during follow-up
  w3.Tr <- sum(data1$Tr == 1 & data1$M1 == 1 & data1$M2 == 1) / denom.Tr    # weight for outcomes imputed with extreme values

  # consistency check
  stopifnot(abs(w1.Tr+w2.Tr+w3.Tr-1) < .00001)

  m1.Tr <- mean(data1$Y[data1$Tr == 1 & data1$M1 == 0])
  m2.Tr <- mean(data1$Y[treated.M2a & data1$M2 == 0])

  # consistency check
  if(length(data1$Y[treated.M2a & data1$M2 == 0]) == 0) m2.Tr <- 0

  # lower bound
  temp[1] <- w1.Tr*m1.Tr + w2.Tr*m2.Tr + w3.Tr*0

  # upper bound
  temp[2] <- w1.Tr*m1.Tr + w2.Tr*m2.Tr + w3.Tr*1

  # calculate less naive estimator that uses randomization and follow-up outcomes but ignores 2nd stage missingness
  temp[6] <- w1.Tr/(w1.Tr+w2.Tr)*m1.Tr + w2.Tr/(w1.Tr+w2.Tr)*m2.Tr

  # calculate Manski bounds for naive estimator
  temp[10] <- m1.Tr*w1.Tr
  temp[11] <- m1.Tr*w1.Tr + (w2.Tr+w3.Tr)


  return(invisible(temp))
  }
