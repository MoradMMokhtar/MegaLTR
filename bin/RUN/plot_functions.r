plot_LTR.chart<-function(data,valuename,groupname,xlabel,ylabel,programname)
{
  point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)
  colnames(data)
  p1 <- ggplot(data=data, aes_string(x=valuename, group=groupname, fill=groupname)) +
    geom_density(adjust=1.5, alpha=.4) + guides(fill=guide_legend(title=ylabel))+
    labs(caption = paste(programname,format(Sys.Date(), "%d %B %Y")),x=xlabel)+
    scale_x_continuous(labels = point,breaks = scales::pretty_breaks(n = 10))+
    scale_y_continuous(labels = point,breaks = scales::pretty_breaks(n = 20))
  p1
}

plot_LTR.boxplot<-function(data,valuename,groupname,xlabel,ylabel,programname)
{
  point <- format_format(big.mark = " ", decimal.mark = ".", scientific = FALSE)
  colnames(data)
  p1 <- ggplot(data=data, aes_string(y=valuename,x=groupname, group=groupname, fill=groupname)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    theme(
      legend.position="none", plot.title = element_text(size=5),
      axis.text.x = element_text(angle = 10, vjust = 1, hjust=1,size = 6),
      axis.title=element_text(size=10,face="bold")
    )+
    labs(caption = paste(programname,format(Sys.Date(), "%d %B %Y")),y=xlabel,x=ylabel)
  p1
}
plot_LTR.save<-function(pname,pplot)
{
  ggsave(paste(pname,".png",sep=""), 
         plot = pplot, 
         device = "png", width = 10, height = 5,
         dpi = 320)
}
