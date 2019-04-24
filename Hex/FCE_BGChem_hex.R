library(ggplot2)
library(ggimage)
library(hexSticker)
library(magick)

fill="white";#"#000000"
color="black";#"#ffffff"
subplot="D:/hex_sticker/fce_BGC/FCE_bgchem_icon.jpg"
size=0.8
s_x=1
s_y=0.95
s_width=0.6
t_x=1
t_y=0.5
text.val="Biogeochemisty\nWorking Group"
t_color="black"
font.val="serif"
t_size=12

hexd <- data.frame(x = 1 + c(rep(-sqrt(3)/2, 2), 0, rep(sqrt(3)/2,2), 0), y = 1 + c(0.5, -0.5, -1, -0.5, 0.5, 1))
hexd <- rbind(hexd, hexd[1, ])
d <- data.frame(x = s_x, y = s_y, image = subplot)
d.txt <- data.frame(x = t_x, y = t_y, label = text.val)

hex.stick=ggplot()+
  geom_polygon(aes_(x = ~x, y = ~y), data = hexd, size = size,fill = fill, color = color)+
  geom_image(aes_(x = ~x, y = ~y, image = ~image),d, size = s_width)+ 
  theme_sticker(size)+
  geom_text(aes_(x = 1, y = 1.6, label = "Biogeochemistry"), size = t_size, color = t_color, family = font.val )+
    geom_text(aes_(x = 1, y = 0.4, label = "Working Group"), size = t_size, color = t_color, family = font.val )+
  geom_text(aes_(x = 1, y = 0.25, label = "FCE LTER"), size = 9, color = rgb(0/255,52/255,102/255), family = font.val,fontface="bold")+  
  geom_text(aes_(x = 0.8, y = 1.4, label = "C"), size = 12, color = adjustcolor("black",0.5), family = font.val,fontface="bold")+
    geom_text(aes_(x = 0.9, y = 0.8, label = "N"), size = 12, color = adjustcolor("black",0.5), family = font.val,fontface="bold")+
    geom_text(aes_(x = 1.4, y = 0.9, label = "P"), size = 12, color = adjustcolor("black",0.5), family = font.val,fontface="bold")
hex.stick


ggsave(hex.stick, width = 43.9, height = 50.8, 
       filename = "D:/hex_sticker/fce_BGC/FCE_BGChem.png", 
       bg = "transparent", units = "mm", dpi=300)

