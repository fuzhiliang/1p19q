#======loh_analysis.oncotarget_v2.bak.pl===========

1)考虑位点间的距离，即window内取占比多的状态;
2)考虑1q、19p的位点，并上图;
3)结合p和q定阈值，即1q/1q <0.7 or 19p/19q <0.7 为阳性，目的是过滤掉整体染色体都缺失的样本;


#======loh_analysis.oncotarget_v2.cnv.bak.pl===========

1)添加了每个臂的cnv值

注:

脚本里是基于GBM panel进行LOH 和CNV 分析，如果进行其他胶质瘤panel,需修改panel文件及baselines文件
