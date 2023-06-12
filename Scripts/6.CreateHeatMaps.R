######################
### Author: Shawn Dove
######################

# This script creates a heat map of trend reliability, a plot of reliability vs weight 
# in the LPI, and a plot of the number of data points in the LPI for each year.


# load packages ----

library(ggplot2)
library(patchwork)
library(ggrepel)


# set directory paths ----

rd_name <- "RDataFiles" # name of directory where .RData files are stored
pd_name <- "Plots" # name of directory where .RData files are stored


# load files ----

realms.results.df <- readRDS(paste(rd_name, "/realms_results_df.RData", sep=""))
terr.results.df <- readRDS(paste(rd_name, "/terr_results_df.RData", sep=""))
fw.results.df <- readRDS(paste(rd_name, "/fw_results_df.RData", sep=""))
marine.results.df <- readRDS(paste(rd_name, "/marine_results_df.RData", sep=""))

# terrestrial realm heat map ----

p1 <- ggplot(data = terr.results.df, 
       aes(x=Taxon, y=Realm, fill=SampPercent100))+
  geom_tile(show.legend=FALSE)+
  geom_text(aes(label=paste(round(SampPercent, 0), "%", sep="")), 
            data=terr.results.df)+
  scale_fill_gradient(low="tan3", high="white")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  ggtitle(label="Terrestrial")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA))


# freshwater heat map ----

p2 <- ggplot(data = FW.results.df, 
       aes(x=Taxon, y=Realm, fill=SampPercent100))+
  geom_tile(show.legend=FALSE)+
  geom_text(aes(label=paste(round(SampPercent, 0), "%", sep="")), 
            data=FW.results.df)+
  scale_fill_gradient(low="turquoise3", high="white")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  ggtitle(label="Freshwater")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA))


# marine heat map ----

p3 <- ggplot(data = Marine.results.df, 
       aes(x=Taxon, y=Realm, fill=SampPercent100))+
  geom_tile(show.legend=FALSE)+
  geom_text(aes(label=paste(round(SampPercent, 0), "%", sep="")), 
            data=Marine.results.df[!is.na(Marine.results.df$SampPercent),])+
  scale_fill_gradient(low="dodgerblue2", high="white")+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))+
  ggtitle(label="Marine")+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA))


# define how the heatmaps are layed out when combined

layout <- "
AABBB
AABBB
#CCC#
#CCC#
"

# combine heatmaps into one plot

p_all <- p1 + p2 + p3 +
  plot_layout(design=layout)

p_all


# save plot

ggsave(filename=paste(pd_name, "/heatmap_update.tiff", sep=""),
       plot = last_plot(),
       device = "tiff",
       width = 10000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")


# plot of trend reliability vs relative weight in LPI ----

ggplot(data = realms.results.df[!is.na(realms.results.df$SampPercent100) & realms.results.df$SampPercent <= 100,], 
       aes(x=RelativeWeight, y=SampPercent100))+
  stat_function(fun=function(x) 100-(1/x))+
  geom_point(size=4)+
  geom_label_repel(aes(label=ifelse((100-SampPercent100)*RelativeWeight>1,paste(System, Realm, Taxon2, sep=" "), "")), 
            size = 5, hjust = 0.5, box.padding=0.25, label.padding=0.35, point.padding=0.3, nudge_x=0.002, nudge_y=3, min.segment.length=0.2)+
  labs(x="Relative weight in LPI", y="Trend reliability (% required populations)")+
  scale_x_continuous()+
  scale_y_continuous(limits=c(0,100))+
  annotate("text", x=0.035, y=65, label="x * (100 - y) = 1", size=6)+
  theme_bw()+
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))


# save plot

ggsave(filename=paste(pd_name, "/reliability_vs_weight_update2.tiff", sep=""),
       plot = last_plot(),
       device = "tiff",
       width = 12000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")

