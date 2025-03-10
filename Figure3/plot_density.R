library(ggplot2)
library(ggnewscale)
library(R.matlab)
library(coin)
library(ggExtra)


# 2Seq, Memory with marginal
#data4 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Ocean13_scope40_peak4_Ent_Mem_pos.mat')
data4 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Groot10_scope40_peak4_Ent_Mem_pos.mat')

df_S1 <- data.frame(x = data4$SPK.S1.pos[,1], y = data4$SPK.S1.pos[,2])
df_S2 <- data.frame(x = data4$SPK.S2.pos[,1], y = data4$SPK.S2.pos[,2])
test_c <- unlist(data4$c.vec)
test_x <- data4$x.vec[,1]
test_y <- data4$y.vec[,1]
df_vec = data.frame(x = test_x, y = test_y, c = test_c)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 3.1,   # Create new plot window
        height = 2.45,
        noRStudioGD = TRUE)
p <- ggplot(data=df_vec,aes(x=x, y=y, color=c)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = c("#6395BC","#CC5F5A")) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) + 
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-22.5) + 
  ylim(-17.5,-5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 2.48,   # Create new plot window
        height = 1.96,
        noRStudioGD = TRUE)
ggplot() + 
  geom_hex(bins=22, data=df_S1, aes(x=x,y=y), alpha = 0.5) +
  scale_fill_gradient(low="snow",high="#6395BC",name="Mem1",na.value="#6395BC",limits=c(0,17)) + 
  new_scale_fill() +
  geom_hex(bins=22, data=df_S2, aes(x=x,y=y), alpha = 0.5) + 
  scale_fill_gradient(low="snow",high="#CC5F5A",name="Mem2",na.value="#CC5F5A",limits=c(0,17)) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) + 
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-22.5) + 
  ylim(-17.5,-5)

# 3Seq, Memory with marginal
#data5 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Ocean13_scope40_peak4_3Seq_Ent_Mem_pos.mat')
data5 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Groot11_scope40_peak4_3Seq_Ent_Mem_pos.mat')

df_S1_3q <- data.frame(x = data5$SPK.S1.pos[,1], y = data5$SPK.S1.pos[,2])
df_S2_3q <- data.frame(x = data5$SPK.S2.pos[,1], y = data5$SPK.S2.pos[,2])
df_S3_3q <- data.frame(x = data5$SPK.S3.pos[,1], y = data5$SPK.S3.pos[,2])
test_c <- unlist(data5$c.vec)
test_x <- data5$x.vec[,1]
test_y <- data5$y.vec[,1]
df_vec = data.frame(x = test_x, y = test_y, c = test_c)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 3.1,   # Create new plot window
        height = 2.45,
        noRStudioGD = TRUE)
p <- ggplot(data=df_vec,aes(x=x, y=y, color=c)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = c("#6395BC","#CC5F5A","#8FBC8F")) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) + 
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-22.5) + 
  ylim(-17.5,-5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 2.48,   # Create new plot window
        height = 1.96,
        noRStudioGD = TRUE)
ggplot() + 
  geom_hex(bins=22, data=df_S1_3q, aes(x=x,y=y), alpha = 0.5) +
  scale_fill_gradient(low="snow",high="#6395BC",name="Mem1",na.value="#6395BC",limits=c(0,17)) + 
  new_scale_fill() +
  geom_hex(bins=22, data=df_S2_3q, aes(x=x,y=y), alpha = 0.5) + 
  scale_fill_gradient(low="snow",high="#CC5F5A",name="Mem2",na.value="#CC5F5A",limits=c(0,17)) + 
  new_scale_fill() +
  geom_hex(bins=22, data=df_S3_3q, aes(x=x,y=y), alpha = 0.5) + 
  scale_fill_gradient(low="snow",high="#8FBC8F",name="Mem3",na.value="#8FBC8F",limits=c(0,17)) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) + 
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-20) + 
  ylim(-17.5,-5)


# 2/3Seq, Entry with marginal
#data6 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Ocean13_scope40_peak4_Ent_Mem_pos.mat')
#data6 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Groot10_scope40_peak4_Ent_Mem_pos.mat')
#data6 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Ocean13_scope40_peak4_3Seq_Ent_Mem_pos.mat')
data6 <- readMat('G:\\PaperPrep\\PostNatNeuroRevision\\SPK_TH_FW_array_dist\\Groot11_scope40_peak4_3Seq_Ent_Mem_pos.mat')

df_TH <- data.frame(x = data6$TH.Ent.pos[,1], y = data6$TH.Ent.pos[,2])
df_SPK <- data.frame(x = data6$SPK.Ent.pos[,1], y = data6$SPK.Ent.pos[,2])
x_vec = c(data6$TH.Ent.pos[,1],data6$SPK.Ent.pos[,1])
y_vec = c(data6$TH.Ent.pos[,2],data6$SPK.Ent.pos[,2])
c_vec = c(replicate(1000,"TH"),replicate(1000,"SPK"))
df_vec = data.frame(x = x_vec, y = y_vec, c = c_vec)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 3.1,   # Create new plot window
        height = 2.45,
        noRStudioGD = TRUE)
p <- ggplot(data=df_vec,aes(x=x, y=y, color=c)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = c("#909090","#F3B664")) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) + 
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-22.5) + 
  ylim(-17.5,-5)
ggMarginal(p, groupColour = TRUE, groupFill = TRUE)

plot.new()            # Create empty plot in RStudios' default window
dev.new(width = 2.48,   # Create new plot window
        height = 1.96,
        noRStudioGD = TRUE)
ggplot() + 
  geom_hex(bins=20, data=df_TH, aes(x=x,y=y), alpha = 0.5) +
  scale_fill_gradient(low="snow",high="#F3B664",name="TH",na.value="#F3B664",limits=c(0,17)) +
  new_scale_fill() +
  geom_hex(bins=20, data=df_SPK, aes(x=x,y=y), alpha = 0.5) + 
  scale_fill_gradient(low="snow",high="#909090",name="SPK",na.value="#909090",limits=c(0,17)) + 
  theme(text = element_text(size = 8),legend.key.size = unit(0.3, 'cm')) +
  # Ocean
  #xlim(-46.5,-35) + 
  #ylim(3,16.5)
  # Groot
  xlim(-37.5,-22.5) + 
  ylim(-17.5,-5)

