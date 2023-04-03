#wczytanie bibliotek
library(readr)
library(dplyr)
library(zoo)
library(ggplot2)
library(rpart)
library(rpart.plot)
library(MASS)
############################### PREPROCESSING ##################################
#import danych
dane <- read.delim("257440.txt", skip = 6, header = FALSE, sep = "")
View(dane)
names(dane) <- c("time","isotope")

#sprawdzenie struktury danych
str(dane) 
count(dane) 
head(dane) 
tail(dane)

#wykonanie wykresów
boxplot(dane$isotope, main="Wykres ramka-wąsy")
plot(dane$isotope, type = "l", main = "Wykres liniowy", xlab = "Indeks", ylab="Stosunek izotopu tlenu")
hist(dane$isotope, breaks=12, main = "Histogram", xlab="Stosunek izotopu tlenu", ylab="Częstość")
summary(dane$isotope)

#sprawdzenie wartości ekstremalnych
srednia <- mean(dane$isotope, na.rm = TRUE)
odchylenie <- sd(dane$isotope, na.rm = TRUE)
dolny_zakres = srednia - 3*odchylenie
gorny_zakres = srednia + 3*odchylenie
dolny_zakres
gorny_zakres

#wykres autokorelacji i autokorelacji cząstkowej
acf(dane$isotope, na.action = na.pass, main="Wykres autokorelacji")
pacf(dane$isotope, na.action = na.pass, main="Wykres autokorelacji cząstkowej")

#dekompozycja szeregu
dane_ts <- ts(dane$isotope, start = 3000, frequency = 100)
dane_dc <- decompose(dane_ts)
plot(dane_dc)

#dodanie zmiennych odzwierciedlających trend liniowy, paraboliczny, hiperboliczny i logarytmiczny
dane_2 <- dane %>%
  mutate(t = c(866:1), t2 = t^2, t3 = t^3, t_log = log(t), 
         sezonowosc = dane_dc$seasonal, trend = dane_dc$trend)
View(dane_2)

#dodanie zmiennej opóźnionej o 11 i o 22,23,24,25
dane_2 <- dane_2 %>%
  mutate(lag_11 = lag((isotope),11))
dane_2 <- dane_2 %>%
  mutate(lag_22 = lag((isotope),22))
dane_2 <- dane_2 %>%
  mutate(lag_23 = lag((isotope),23))
dane_2 <- dane_2 %>%
  mutate(lag_24 = lag((isotope),24))
dane_2 <- dane_2 %>%
  mutate(lag_25 = lag((isotope),25))

#korelacja
cor(dane_2$lag_11, dane_2$isotope, use="complete.obs")
cor(dane_2$lag_22, dane_2$isotope, use="complete.obs") 
cor(dane_2$lag_23, dane_2$isotope, use="complete.obs")
cor(dane_2$lag_24, dane_2$isotope, use="complete.obs")
cor(dane_2$lag_25, dane_2$isotope, use="complete.obs")

library(corrplot)
x <- cor(dane_2, method = "pearson", use="complete.obs")
corrplot(x)
x <- cor(dane_2, method = "spearman", use="complete.obs")
corrplot(x)


######################### DRZEWO REGRESYJNE CART ###############################
library(rpart)
library(rpart.plot)

#podział na zbiór uczący i testowy
uczacy <- dane_2[-c(863, 864, 865, 866),]
testowy <- dane_2[c(863, 864, 865, 866),]

#model
set.seed(10)
rt_dane <- rpart(isotope~t+t2+t3+t_log+lag_11+lag_22+lag_23+lag_24+lag_25, 
                 data = uczacy, control = rpart.control(cp=0.00001))
printcp(rt_dane)

#wykres drzewa przed przycinaniem
plot(rt_dane)
text(rt_dane, cex = 0.9, xpd = TRUE)

#przycinanie drzewa i jego wizualizacja
plotcp(rt_dane)
rt_dane_pr<-prune(rt_dane, cp=0.006)
plot(rt_dane_pr)
text(rt_dane_pr, cex=0.9, xpd=TRUE)

rpart.plot(rt_dane_pr)

#wygenerowany model i reszty, ocena jakości modelu
rt_dane_pr_mod <-predict(rt_dane_pr)
R2 <- cor(uczacy$isotope, rt_dane_pr_mod)^2
rt_dane_pr_res <- uczacy$isotope - rt_dane_pr_mod
blad = mean(uczacy$isotope - rt_dane_pr_mod)
plot(uczacy$isotope, type="l", col="black", main = "Dopasowanie modelu do danych")
lines(rt_dane_pr_mod, type="l", col="red")
#printcp(rt_dane)

#ocena jakości reszt
hist(rt_dane_pr_res, main = "Histogram reszt") 
qqnorm(rt_dane_pr_res, pch = 1, frame = FALSE, main = "Wykres normalności")
qqline(rt_dane_pr_res, col = "steelblue", lwd = 2)
acf(rt_dane_pr_res, main = "Wykres autokorelacji reszt") 
pacf(rt_dane_pr_res, main = "Wykres autokorelacji cząstkowej reszt")
plot(rt_dane_pr_res, main = "Wykres rozrzutu reszt")


#przewidywanie na zbiór testowy
library(MLmetrics)
rt_dane_pr_pred <- predict(rt_dane_pr, testowy[,-2])

#MAPE, THEIL
MAPE(rt_dane_pr_pred, testowy$isotope)
library(REAT)
theil(rt_dane_pr_pred, na.rm = TRUE)



####################### MODEL REGRESJI WIELORAKIEJ #############################

#podział na zbiór uczący i testowy
uczacy <- dane_2[-c(863, 864, 865, 866),]
testowy <- dane_2[c(863, 864, 865, 866),]

#wykonanie modelu
m <- lm(isotope~t+t2+t3+t_log+lag_11+lag_22+lag_23+lag_24+lag_25, data = uczacy)

#ocena jakości modelu
summary(m)
AIC(m)
BIC(m)

#dopasowanie modelu do danych
plot(uczacy$isotope, type="l", col="black", main = "Dopasowanie modelu do danych")
lines(m$fitted.values, type="l", col="red")

#ocena jakości reszt
hist(m$residuals, main = "Histogram reszt") 
qqnorm(m$residuals, pch = 1, frame = FALSE, main = "Wykres normalności")
qqline(m$residuals, col = "steelblue", lwd = 2)
acf(m$residuals, main = "Wykres autokorelacji reszt") 
pacf(m$residuals, main = "Wykres autokorelacji cząstkowej reszt")
plot(m$residuals~m$fitted.values, main = "Wykres rozrzutu reszt")

#ocena jakości prognoz
library(MLmetrics)

pred_m<-predict(m, newdata = testowy)
MAPE(pred_m, testowy$isotope)
theil(pred_m, na.rm = TRUE)



#nieznaczne poprawienie modelu regresji wielorakiej
m <- lm(isotope~t2+t3+lag_11+lag_22+lag_25, uczacy)
summary(m) #wszystkie zmienne istotne statystycznie
AIC(m)
BIC(m)

plot(uczacy$isotope, type="l", col="red")
lines(m$fitted.values, type="l", col="blue")

hist(m$residuals) 
qqnorm(m$residuals, pch = 1, frame = FALSE)
qqline(m$residuals, col = "steelblue", lwd = 2)
acf(m$residuals) 
pacf(m$residuals)
plot(m$residuals~m$fitted.values) 

pred_m<-predict(m, newdata = testowy)
MAPE(pred_m, testowy$isotope) #poprawa mape
theil(pred_m, na.rm = TRUE)

