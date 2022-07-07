###XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX###
### Fit models to compare habitat suitability and functionality ###
###XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX###

# * Data load ============
library(raster)
load("data/datastack_sum.rda")
datastack

# ** Coarse grain & convert to dataframe #############
coarse_stack <- aggregate(datastack, fact=5, fun=mean, na.rm=T)
coarse_stack

plot(coarse_stack[[6]])
coarsedf <- as.data.frame(coarse_stack,xy=TRUE)
head(coarsedf)
coarsedf <- na.omit(coarsedf)

tmpdf <- coarsedf

# * Fit monotonic gam =====
library(scam)

names(tmpdf)

wcol <- c(4:ncol(tmpdf))
wcol <- wcol[!(wcol %in% grep("use", names(tmpdf)))]
wcol <- wcol[!(wcol %in% grep("season", names(tmpdf)))]
wcol <- wcol[!(wcol %in% grep("winter", names(tmpdf)))]
(wvar <- names(tmpdf)[wcol])

head(table_fit <- data.frame(var=wvar, aic=NA, rsq=NA, 
                             dev_expl=NA, SSR=NA))

i=1
for (i in c(1:length(wvar))){ 
  ffm <- as.formula(paste0('use ~ s(', wvar[i], ',k=5, bs="mpi")'))
  gam_mod <- scam(formula=ffm, data=tmpdf)
  
  table_fit[i, c("aic")] <- AIC(gam_mod)
  table_fit[i, c("rsq")] <- summary(gam_mod)$r.sq
  table_fit[i, c("dev_expl")] <- summary(gam_mod)$dev.expl
  
  pred_vals <- predict(gam_mod, newdata=tmpdf, type="response")
  table_fit$SSR[i] <- sum((tmpdf$use - pred_vals)^2)
}

table_fit

table_fit[order(table_fit$aic),]
table_fit[order(-table_fit$rsq),]
table_fit[order(table_fit$SSR),]


# ** prep table ##########
table_fit$dev_expl <- round(table_fit$dev_expl,2) 
table_fit$daic <- table_fit$aic-min(table_fit$aic) 

table_fit$aic_weight <- exp(-0.5*table_fit$daic)
table_fit$aic_weight <- table_fit$aic_weight/sum(table_fit$aic_weight) 
table_fit$daic <- round(table_fit$daic,0) 
table_fit$aic_weight <- round(table_fit$aic_weight,2) 

table_fit[order(table_fit$aic),]
table_fit[order(table_fit$daic),-c(2,3,5)]

ffm <- as.formula(paste0('use ~ s(func_params01_3,k=5, bs="mpi")'))
gam_mod <- scam(formula=ffm, data=tmpdf)
summary(gam_mod)
plot(gam_mod)

test <- residuals(gam_mod)
hist(test, breaks=100)

ffm <- as.formula(paste0('use ~ s(quality,k=5, bs="mpi")'))
gam_mod2 <- scam(formula=ffm, data=tmpdf)
plot(gam_mod2)

test <- residuals(gam_mod2)
hist(test, breaks=100)

# * Predict models =====
library(ggplot2)
library(reshape2)

mod1 <- scam(use~s(func_params01_3, k=5, bs="mpi"),data=tmpdf)
mod1 <- scam(use~s(func_params01_3, k=5, bs="mpi"),data=tmpdf)
plot(mod1)

mod2 <- scam(use~s(quality, k=5, bs="mpi"),data=tmpdf)
mod2 <- scam(use~s(quality, k=10, bs="mpi"),data=tmpdf)
plot(mod2)


tmpdf$stdz_func <- tmpdf$func_params01_3/quantile(tmpdf$func_params01_3, probs=0.99) 
tmpdf$stdz_func <- ifelse(tmpdf$stdz_func>1, 1, tmpdf$stdz_func) 
tmpdf$stdz_quality <- tmpdf$quality/quantile(tmpdf$quality, probs=0.99) 
tmpdf$stdz_quality <- ifelse(tmpdf$stdz_quality>1, 1, tmpdf$stdz_quality) 

tmpdf$stdz_func2 <- round(tmpdf$stdz_func, 1) 
tmpdf$stdz_quality2 <- round(tmpdf$stdz_quality, 1) 


plt <- ggplot(tmpdf, aes(x=factor(stdz_quality2), y=use)) + 
  geom_boxplot(col="#74D055FF", fill=alpha("#74D055FF", 0.5), outlier.size=0.5, outlier.alpha=0.25) + 
  xlab("scaled suitability") + ylab("reindeer distribution") + ylim(0,0.75) +
  theme_classic()
plt

plt2 <- ggplot(tmpdf, aes(x=factor(stdz_func2), y=use)) + 
  geom_boxplot(col="#1F968BFF", fill=alpha("#1F968BFF", 0.5), outlier.size=0.5, outlier.alpha=0.25) + 
  xlab("scaled functionality") + ylab("reindeer distribution") + ylim(0,0.75) +
  theme_classic()
plt2

require(gridExtra)
grid.arrange(plt, plt2, nrow=2)

tmp <- data.frame(stdz_quality2=c(0:10)/10)
tmp$quality <- tmp$stdz_quality2*quantile(tmpdf$quality, probs=0.99)
tmp$mn <- predict(mod2, newdata=tmp, se.fit=T)[[1]]
tmp$se <- predict(mod2, newdata=tmp, se.fit=T)[[2]]

tmp$upper <- tmp$mn + 2*tmp$se
tmp$lower <- tmp$mn - 2*tmp$se
tmp$mn <- ifelse(tmp$mn<0, 0, tmp$mn)
tmp$upper <- ifelse(tmp$upper<0, 0, tmp$upper)
tmp$lower <- ifelse(tmp$lower<0, 0, tmp$lower)

plt <- plt + geom_line(data=tmp, aes(x=factor(stdz_quality2), y=mn, group=1), col="#74D055FF", lwd=1.5) +
  geom_ribbon(data=tmp, aes(x=factor(stdz_quality2), y=mn, ymin = lower, ymax = upper, group=1), col="#74D055FF", alpha=0.25)
plt

tmp2 <- data.frame(stdz_func2=c(0:10)/10)
tmp2$func_params01_3 <- tmp2$stdz_func2*quantile(tmpdf$func_params01_3, probs=0.99)
tmp2$mn <- predict(mod1, newdata=tmp2, se.fit=T)[[1]]
tmp2$se <- predict(mod1, newdata=tmp2, se.fit=T)[[2]]

tmp2$upper <- tmp2$mn + 2*tmp2$se
tmp2$lower <- tmp2$mn - 2*tmp2$se
tmp2$mn <- ifelse(tmp2$mn<0, 0, tmp2$mn)
tmp2$upper <- ifelse(tmp2$upper<0, 0, tmp2$upper)
tmp2$lower <- ifelse(tmp2$lower<0, 0, tmp2$lower)

plt2 <- plt2 + geom_line(data=tmp2, aes(x=factor(stdz_func2), y=mn, group=1), col="#1F968BFF", lwd=1.5) +
  geom_ribbon(data=tmp2, aes(x=factor(stdz_func2), y=mn, ymin = lower, ymax = upper, group=1), col="#1F968BFF", alpha=0.25)
plt2

plt <- grid.arrange(plt, plt2, nrow=2)
plt

tmpdf2 <- rbind(data.frame(use=tmpdf$use, variable="functionality", value=tmpdf$stdz_func2), 
                data.frame(use=tmpdf$use, variable="suitability", value=tmpdf$stdz_quality2))

plt <- ggplot(tmpdf2, aes(x=factor(value), y=use)) + 
  geom_boxplot(aes(col=factor(variable, levels = c("suitability", "functionality")), fill=factor(variable, levels = c("suitability", "functionality"))), outlier.size=0.5, outlier.alpha=0.25) +
  scale_fill_manual(name="predictor", labels=c("suitability", "functionality"), values=alpha(c("#74D055FF", "#1F968BFF"), 0.5)) +
  scale_color_manual(name="predictor", labels=c("suitability", "functionality"), values=c("#74D055FF", "#1F968BFF")) + 
  xlab("scaled suitability/functionality value") + ylab("reindeer distribution") + ylim(0,0.75) +
  theme_classic()
plt

plt <- plt + geom_line(data=tmp, aes(x=factor(stdz_quality2), y=mn, group=1), col="#74D055FF", lwd=1.5) +
  geom_ribbon(data=tmp, aes(x=factor(stdz_quality2), y=mn, ymin = lower, ymax = upper, group=1), col="#74D055FF", alpha=0.25)
plt

plt <- plt + geom_line(data=tmp2, aes(x=factor(stdz_func2), y=mn, group=1), col="#1F968BFF", lwd=1.5) +
  geom_ribbon(data=tmp2, aes(x=factor(stdz_func2), y=mn, ymin = lower, ymax = upper, group=1), col="#1F968BFF", alpha=0.25)
plt

ggsave("figs/shape_summer.png", plt, width = 20, height = 15, units = "cm", dpi="print")
