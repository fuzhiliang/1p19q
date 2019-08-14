library(ggplot2)
data <- read.table("test/sample20190412-B2-A00168_capture20190411-2_I1/sample20190412-B2-A00168_capture20190411-2_I1.snv.chr1_chr19.ratio.txt",header = F)
data1 <- data[,c(2,3,8)]
fish <- data.frame(chr = c("chr1","chr1", "chr19","chr19"),  pos = c(0, 7.2, 43.4,51.4))
pq <- data.frame(chr = c("chr1", "chr19"),  pos = c(125,26.5))
names(data1) <- c("chr","pos","ratio")
df <- data.frame(x = data1$pos,y=data1$ratio,chr=data1$chr)
p <- ggplot(df,aes(x=x/1000000,y=y)) + geom_point(size = 0.9) + facet_wrap(~chr,scales = "free",nrow = 2) + xlab("pos") + ylab("ratio") + labs(title="sample20190412-B2-A00168_capture20190411-2_I1") + theme(plot.title = element_text(hjust = 0.5)) + scale_color_manual(values = "red") + scale_y_continuous(limits = c(0,1))  + geom_vline(aes(xintercept = pos),color = "red", fish) + geom_vline(aes(xintercept = pos) , color = 'red', size = 0.2, linetype=2,pq)+ geom_hline(yintercept = 0.4 , color = 'red', size = 0.2, linetype="solid") + geom_hline(yintercept = 0.6 , color = 'red', size = 0.2, linetype="solid") 
ggsave(p, file = "test/sample20190412-B2-A00168_capture20190411-2_I1/sample20190412-B2-A00168_capture20190411-2_I1.snv.ratio.png", width = 12, height = 8, dpi = 100)
