library(ggplot2)
library(patchwork)
library(ggrepel)

# generate bar charts
p1 <- ggplot(data = realms.results.df[realms.results.df$System=="Terrestrial",], 
       aes(y=Taxon, x=SampPercent100, fill=Taxon))+
  geom_bar(stat='identity', show.legend=FALSE)+
  facet_grid(Realm~System, scales="free")+
  theme(axis.title.y=element_blank(),
        axis.title.x=element_blank())

p2 <- ggplot(data = realms.results.df[realms.results.df$System=="Freshwater",], 
             aes(y=Taxon, x=SampPercent100, fill=Taxon))+
  geom_bar(stat='identity', show.legend=FALSE)+
  facet_grid(Realm~System, scales="free")+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks=element_blank())
  
p3 <- ggplot(data = realms.results.df[realms.results.df$System=="Marine",], 
             aes(y=Taxon, x=SampPercent100, fill=Taxon))+
  geom_bar(stat='identity')+
  facet_grid(Realm~System, scales="free")+
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks=element_blank())

p_all <- p1 + p2 + p3
p_all

p1 <- ggplot(data = Terr.results.df, 
       aes(x=Taxon, y=Realm, fill=SampPercent100))+
  geom_tile(show.legend=FALSE)+
  geom_text(aes(label=paste(round(SampPercent, 0), "%", sep="")), 
            data=Terr.results.df)+
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
        #axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill=NA))

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

layout <- "
AABBB
AABBB
#CCC#
#CCC#
"

p_all <- p1 + p2 + p3 +
  plot_layout(design=layout)#+
  #plot_annotation(title = "% of minimum sample size achieved",
  #                theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))


p_all

ggsave("heatmap.tiff",
       plot = last_plot(),
       device = "tiff",
       width = 10000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")


ggplot(data = realms.results.df[!is.na(realms.results.df$SampPercent100) & realms.results.df$SampPercent <= 100,], 
       aes(x=RelativeWeight, y=SampPercent100))+
  stat_function(fun=function(x) 100-(1/x))+
  geom_point(size=4)+
  geom_label_repel(aes(label=ifelse((100-SampPercent100)*RelativeWeight>1,paste(System, Realm, Taxon2, sep=" "), "")), 
            size = 4, hjust = 0, box.padding=0.25, label.padding=0.35, point.padding=0.3, nudge_x=0.002, nudge_y=3, min.segment.length=0.2)+
  labs(x="Relative weight in LPI", y="Trend reliability (% required populations)")+
  scale_x_continuous()+
  scale_y_continuous(limits=c(0,100))+
  annotate("text", x=0.055, y=77, label="x * (100 - y) = 1", size=6)+
  theme_bw()+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16, hjust = 0.5),
        panel.background = element_rect(fill="white"))

ggsave("reliability_vs_weight.tiff",
       plot = last_plot(),
       device = "tiff",
       width = 12000,
       height = 8000,
       units = "px",
       dpi = 1000,
       compression = "lzw")

#########
#Pearsons correlation coefficient testing
pcdf <- realms.results.df[!is.na(realms.results.df$SampPercent) & !is.na(realms.results.df$RelativeWeight),]
pcdf.ter <- pcdf[pcdf$System=="Terrestrial",]
pcdf.fw <- pcdf[pcdf$System=="Freshwater",]
pcdf.mar <- pcdf[pcdf$System=="Marine",]
summary(cor(pcdf$SampPercent, pcdf$RelativeWeight, method="pearson"))
cor.test(pcdf.ter$SampPercent, pcdf.ter$RelativeWeight)
cor.test(pcdf.fw$SampPercent, pcdf.fw$RelativeWeight)
cor.test(pcdf.mar$SampPercent, pcdf.mar$RelativeWeight)
