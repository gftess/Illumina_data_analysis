#########   ANÁLISE DOS DADOS DO SEQUENCIAMENTO ILLUMINA DA BIBLIOTECA BISON ########

#### Agrupamento da contagem dos sgRNAs de acordo com o gene alvo ####

### O agrupamento foi feito tanto para as contagens da representação inicial (og), 
### para a contagem da biblioteca do DNA plasmidial (p) e para a da biblioteca 
### do RNA lentiviral (l)

sum_count <- tapply(count_gene$count, count_gene$Target_Gene_Symbol, FUN=sum) ## vetor desorganizado
df_sum_count <- data.frame(sum_count) ## transformar em dataframe 
order_sum <- df_sum_count[order(df_sum_count$sum_count),] ## vetor organizado 

#### Ordem da menor contagem para a maior ####

## primeiro ordenar em dataframe
order_og_df <- bison_og[order(bison_og$Library_Representation),]
order_count_p <- count_reads_plasmidial[order(count_reads_plasmidial$count),] 
order_count_l <- count_reads_lentiviral[order(count_reads_lentiviral$count),] 

### criar um vetor só com o valor dos counts ordenados ###
order_og <- order_og_df$Library_Representation
order_p <- order_count_p$count   
order_l <- order_count_l$count

order_l[is.na(order_l)] <- 0  ## substituir os NAs por zero

### Cálculo da média e desvio padrão + minímo e máximo ###
mean()
sd()
min()
max()

### Definir o valor de dois Desvios padrão para formação dos intervalos de 
### interesse 

#### Barplots tanto da contagem de cada sgRNA como do agrupamento dos genes ####

## Plasmidial ## Entidades ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_p, c(-Inf,182.9144 ,1782.234, Inf) ) ]
### o intervalo da média -/+ 2 sd terá a cor "darkred" 
### enquanto as contagens fora do intervalo serão "grey70"

barplot_p <- barplot(order_p, border = NA, width = 0.01,
                     col = cols6, ylab=substitute(paste(italic('Reads Illumina'))) ,
                     xlab = "Entidades da biblioteca Bison", 
                     xlim = c(0,36), ylim = c(0,3100))

abline(h=982.5744, col= "black", lwd = 1.75) ## média
abline(h=91, col= "grey60", lwd = 1.75, lty = "dashed") ## minímo 
abline(h=2970, col= "grey60", lwd = 1.75, lty = "dashed") ## máximo 
box()

op_min <- order_p>182.9144
sum(op_min) # 7
op_max <- order_p < 1782.234
sum(op_max) #110   ######## 95.89474% em vermelho 

## Plasmidial ## Genes ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_sum_p, c(-Inf,2335.08 ,5529.92, Inf) ) ]

barplot_test3 <- barplot(order_sum_p, border = NA, width = 0.01, col = cols6, ylim = c(0,8500),
                         ylab =substitute(paste(italic('Reads Illumina'))),
                         xlab = "Genes alvos da biblioteca Bison" )

abline(h=3932.503, col= "black", lwd = 1.75) ## média
abline(h=1990, col= "grey60", lwd = 1.75, lty = "dashed") ## minímo 
abline(h=8401, col= "grey60", lwd = 1.75, lty = "dashed") ## máximo 
box()

op_min <- order_sum_p1 > 2335.08
sum(op_min) # 7
op_max <- order_sum_p1 < 5529.92
sum(op_max) #22   ######## 95.92697% em vermelho 

## Lentiviral ## Entidades ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_l, c(-Inf,280.005 ,1750.405, Inf) ) ]


barplot_l <- barplot(order_l, ylim=c(0,3000), border = NA, width = 0.01, col = cols,
                     ylab=substitute(paste(italic('Reads Illumina'))) , xlab = "Entidades da biblioteca Bison",
                     xlim = c(0,37))

abline(h=1015.205, col= "black", lwd = 1.75) ## média
abline(h=97, col= "grey70", lwd = 1.75, lty = "dashed") ## minímo 
abline(h=2787, col= "grey70", lwd = 1.75, lty = "dashed") ## máximo 
box()

ol_min <- order_l1 < 280.005
sum(ol_min) # 14
ol_max <- order_l1 > 1750.405
sum(ol_max) #96   ######## 96.139% em vermelho


## Lentiviral ## Genes ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_sum_l, c(-Inf,2592.923 ,5527.183, Inf) ) ]


barplot_sl <- barplot(order_sum_l3, border = NA, width = 0.01, col = cols3, ylim = c(0,8000),
                      ylab =substitute(paste(italic('Reads Illumina'))),
                      xlab = "Genes alvos da biblioteca Bison" )

abline(h=4060.053, col= "black", lwd = 1.75) ## média
abline(h=2469, col= "grey70", lwd = 1.75, lty = "dashed") ## minímo 
abline(h=7591, col= "grey70", lwd = 1.75, lty = "dashed") ## máximo 
box()

ol_min <- order_sum_l3 > 2592.923
sum(ol_min) # 8
ol_max <- order_sum_l3 < 5527.183
sum(ol_max) #23   ######## 95.65% em vermelho

## Bison ## Entidades ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_og, c(-Inf,119.97 ,489.69, Inf) ) ]

barplot_og <- barplot(order_og, border = NA, width = 0.01, col = cols6,
                      ylab ="Representação Inicial",
                      xlab = "Entidades da biblioteca Bison", ylim = c(0, 700) )

abline(h=304.83, col= "black", lwd = 1.75) ## média
abline(h=31.35, col= "grey70", lwd = 1.75, lty = "dashed") ## minímo ## 
abline(h=694, col= "grey70", lwd = 1.75, lty = "dashed") ## máximo ##
box()

og_min <- order_og > 119.97
sum(og_min) # 32
og_max <- order_og < 489.69
sum(og_max) #89   ######## 95.75736% em vermelho


## Bison ## Genes ##

cols6   <- c("grey70", "darkred", "grey70") [findInterval(order_sum_og, c(-Inf, 835.203 ,1603.403, Inf) ) ]

barplot_sog <- barplot(order_sum_og, border = NA, width = 0.01, col = cols6, ylim = c(0,2000),
                       ylab ="Representação Inicial",
                       xlab = "Genes alvos da biblioteca Bison" )

abline(h=1219.303, col= "black", lwd = 1.75) ## média
abline(h=694, col= "grey70", lwd = 1.75, lty = "dashed") ## minímo ## 
abline(h=1916.5, col= "grey70", lwd = 1.75, lty = "dashed") ## máximo ##
box()



og_min <- order_sum_og > 835.203
sum(og_min) # 8
og_max <- order_sum_og < 1603.403
sum(og_max) #19   ######## 96.21318% em vermelho



#### Piecharts tanto da contagem de cada sgRNA como do agrupamento dos genes #####
cols9 <- hcl.colors(4, palette = "Blues", rev = FALSE)

## Plasmidial ## Entidades ##

## porcentagem dos intervalos ### 
# < 2 SD = 0 - 182.9144
7/2850 * 100 ### 0.245614 %
# Média - 2 SD = 182.9144 - 982.74 
1565/2850 * 100 ### 54.91228%
# Média + 2 SD = 982.74 - 1782.234
1168/2850 * 100 ### 40.98246%
# > 2 SD = > 1782.234
110/2850 * 100 ### 3.859649%

percen_p <- c(0.245614, 54.91228, 40.98246,3.859649)
pie(percen_p, labels = c("0.25%", "54.91%", "49.98%", "3.86%"), cex = 0.7 , col = cols9)
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = cols9, cex = 0.75, title = "Desvio Padrão", title.font = 1)

## Plasmidial ## Genes ##

## porcentagem dos intervalos ###
# < 2 SD =  0 - 2335.08
7/712 * 100 ### 0.9831461 %
# Média - 2 SD = 2335.08 - 3932.50 
373/712 * 100### 52.38764%
# Média + 2 SD = 3932.50 - 5529.92
310/712 * 100 ### 43.53933%
# > 2 SD =  > 5529.92
22/712 * 100 ### 3.089888%

percen_p <- c(0.9831461, 52.38764, 43.53933,3.089888)
pie(percen_p, labels = c("0.98%", "52.39%", "43.54%", "3.09%"), cex = 0.7 ,
    col = palheta(4), main = "Genes (plasmidial)")
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = palheta(4), cex = 0.75, title = "Desvio Padrão", title.font = 1)


## Lentiviral ## Entidades ##

## porcentagem dos intervalos ### 
# < 2 SD = 0 - 280.005
14/2849 * 100 ### 0.4914005 %
# Média - 2 SD = 280.005 - 1015.205 
1516/2849 * 100### 53.21165 %
# Média + 2 SD = 1015.205 - 1750.405
1223/2849 * 100 ### 42.92734 %
# > 2 SD = > 1750.405
96/2849 * 100 ### 3.369603 %

percen_p <- c(0.4914005, 53.21165, 42.92734, 3.369603)
pie(percen_p, labels = c("0.49%", "53.21%", "42.93%", "3.37%"), cex = 0.7 , col = cols9)
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = cols9, cex = 0.75, title = "Desvio Padrão", title.font = 1)

## Lentiviral ## Genes ##

## porcentagem dos intervalos ### 
# < 2 SD = 0 - 2592.923
8/712 * 100 ### 1.123596 %
# Média - 2 SD = 2592.923 - 4060.05 
370/712 * 100### 51.96629 %
# Média + 2 SD = 4060.05  - 5527.183
311/712 * 100 ### 43.67978 %
# > 2 SD = > 5527.183
23/712 * 100 ### 3.230337 %

percen_p <- c(1.123596,51.96629 , 43.67978, 3.230337)
pie(percen_p, labels = c("1.12%", "51.97%", "43.68%", "3.23%"), cex = 0.7 , col = palheta(4),
    main = "Genes (lentiviral)")
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = palheta(4), cex = 0.75, title = "Desvio Padrão", title.font = 1)


## Bison ## entidades ##

## porcentagem dos intervalos ### 
# < 2 SD = 0 - 119.97
32/2852 * 100 ### 1.12202 %
# Média - 2 SD = 119.97 - 304.83
1486/2852 * 100### 52.10379 %
# Média + 2 SD = 304.83  - 489.69
1245/2852 * 100 ### 43.65358 %
# > 2 SD = > 489.69
89/2852 * 100 ### 3.120617 %

percen_p <- c(1.12202, 52.10379 ,43.65358 , 3.120617)
pie(percen_p, labels = c("1.12%", "52.10%", "43.66%", "3.12%"), cex = 0.7 , col = cols9)
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = cols9, cex = 0.75, title = "Desvio Padrão", title.font = 1)


## Bison ## Genes ##

## porcentagem dos intervalos ### 
# < 2 SD = 0 - 835.203
8/713 * 100 ### 1.12202 %
# Média - 2 SD = 835.203 - 1219.30
354/713 * 100### 49.64937 %
# Média + 2 SD = 1219.30  - 1603.403
332/713 * 100 ### 46.56381 %
# > 2 SD = > 1603.403
19/713 * 100 ### 2.664797 %

percen_p <- c(1.12202,49.64937 ,46.56381, 2.664797)
pie(percen_p, labels = c("1.12%", "49.65%", "46.56%", "2.67%"), cex = 0.7 , col = palheta(4),
    main = "Genes (inicial)")
legend("bottomleft",legend=c("< 2 sd", "Média - 2 sd", "Média + 2 sd", "> 2 sd"),
       fill = palheta(4), cex = 0.75, title = "Desvio Padrão", title.font = 1)



#### Gráficos de correlação entre a contagem dos genes inicial e entre as bibliotecas #####

# Plasmidial ## Genes ##
## foi feita uma dataframe alinhando os genes por ordem alfabética
## com a contagem das três bibliotecas 

p_plot <- plot(reglim$sum_og, reglim$sum_count_p, xlim = c(0,2500),
               ylim = c(0, 9000), ylab = "Contagem dos genes (plasmidial)",
               xlab = "Contagem dos genes (inicial)", 
               col = "lightsteelblue", pch = 19) ## plot com os pontos


correlation_p <- cor(reglim$sum_og, reglim$sum_count_p, method = 'pearson') # correlação de Pearson
correlation_p
reaglin <- lm(reglim$sum_count_p~reglim$sum_og)
abline(reaglin, col = "red4", lty = "dashed")


## Lentiviral ## Genes ##

l_plot <- plot(reglim$sum_og, reglim$sum_count_l, xlim = c(0,2500),
               ylim = c(0, 9000), ylab = "Contagem dos genes (lentiviral)",
               xlab = "Contagem dos genes (inicial)", 
               col = "steelblue4", pch = 19) ## plot com os pontos


correlation_l <- cor(reglim$sum_og, reglim$sum_count_l, method = 'pearson') # correlação de Pearson
correlation_l
reaglin_l <- lm(reglim$sum_count_l~reglim$sum_og)
abline(reaglin_l, col = "red4", lty = "dashed")
