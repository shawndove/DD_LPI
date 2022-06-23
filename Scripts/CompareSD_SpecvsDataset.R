counter <- 1
temp2 <- list()
for (i in 1:length(temp)) {
  if (length(temp[[i]]) > 2) {
    
    temp2[[counter]] <- temp[[i]]
    counter <- counter + 1
  }
}

temp3 <- list()
for (i in 1:length(temp2)) {
 temp3[[i]] <- mean(temp2[[i]])
}

temp2_sd <- vector()
for (i in 1:length(temp2)) {
  temp2_sd[i] <- sd(temp2[[i]])
}


