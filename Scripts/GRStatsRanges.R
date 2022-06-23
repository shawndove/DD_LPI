meangr <- vector()
sdgr <- vector()
meansd <- vector()
for (i in 1:length(gr.stats.list)) {
  

  if (!is.na(pop.size.vec[i]) & pop.size.vec[i] >= 100) {
    
    meangr[i] <- gr.stats.list[[i]][1]
    sdgr[i] <- gr.stats.list[[i]][2]
    meansd[i] <- gr.stats.list[[i]][3]
    
  } else {
    
    meangr[i] <- NA
    sdgr[i] <- NA
    meansd[i] <- NA
    
  }

  
}
meangr <- unlist(meangr)
sdgr <- unlist(sdgr)
meansd <- unlist(meansd)

mean(meangr, na.rm=TRUE)
mean(sdgr, na.rm=TRUE)
mean(meansd, na.rm=TRUE)
range(meangr[pop.size.vec >20 ], na.rm=TRUE)
range(sdgr[pop.size.vec >20 ], na.rm=TRUE)
range(meansd[pop.size.vec >20 ], na.rm=TRUE)

range(mean.tslength.vec[pop.size.vec >20 ], na.rm=TRUE)
