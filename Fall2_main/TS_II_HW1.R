library(tseries)
library(aTSA)
library(forecast)
library(tidyverse)
library(ggplot2)
library(readr)
library(httpgd)
library(zoo)

setwd("Homework1_TS2")
hgd(port = 8000)

# Read in the training and validation data
# Combine the two current validation data into one
energy <- read_csv("hrl_load_metered.csv")
valid1 <- read_csv("hrl_load_metered_test1.csv")
valid2 <- read_csv("hrl_load_metered_test2.csv")
valid <- bind_rows(valid1, valid2)

# Checking for duplicate hours in the train and valid data due to daylight savings
energy %>%
    filter(duplicated(datetime_beginning_ept) | duplicated(datetime_beginning_ept, fromLast = TRUE))

valid %>%
    filter(duplicated(datetime_beginning_ept) | duplicated(datetime_beginning_ept, fromLast = TRUE))


dl_sav <- energy$datetime_beginning_ept %in% names(table(energy$datetime_beginning_ept)[table(energy$datetime_beginning_ept) > 1])
impute <- energy[which(dl_sav) - 7*34*24 + 1,]

impute <- impute %>% rowwise( ) %>% mutate(date = str_split(datetime_beginning_ept, " ")[[1]][1], hour = as.numeric(str_split(str_split(datetime_beginning_ept, ":")[[1]][1], " ")[[1]][2]) )

impute <- impute %>% group_by(date) %>% summarize(mw = mean(mw))
impute <- impute %>% rowwise() %>% mutate(splitted = str_split(date, "/")) %>% mutate(month = as.numeric(splitted[[1]][1]), day = as.numeric(splitted[[2]][1]), year = as.numeric(splitted[[3]][1]), hour = 2) %>% select(year, month, day, hour, mw)

# Because of daylight savings, some hours are missing and some are duplicated. We're going to impute values for the missing hours and compress duplicated hours into a single average.
temp <- energy %>% rowwise() %>% mutate(date = str_split(datetime_beginning_ept, " ")[[1]][1], hour = as.numeric(str_split(str_split(datetime_beginning_ept, ":")[[1]][1], " ")[[1]][2]) ) %>% mutate(splitted = str_split(date, "/")) %>% mutate(day = as.numeric(splitted[[2]][1]), year = as.numeric(splitted[[3]][1]), month = as.numeric(splitted[[1]][1])) %>% select(year, month, day, hour, mw)

temp <- rbind(temp, impute)
temp <- temp %>% group_by(year, month, day, hour) %>% summarize(mw = mean(mw))
temp <- temp[order(temp$year, temp$month, temp$day, temp$hour),]
energy <- temp

# Create time series data with seasonality length S = 24 (frequency)
# Frequency represents the seasonal length with respect to unit of time expressed in data
energy_ts <- ts(energy$mw, start = c(1, 0), frequency = 24)

# Plot overall time series data
# There appears to be a seasonality component throughout so we should address before modeling ARIMA
autoplot(energy_ts)

# Create add and mul Holt-Winters model that forecasts over the validation period
energy_hw_add <- hw(energy_ts, seasonal = "additive", h = nrow(valid))
energy_hw_mul <- hw(energy_ts, seasonal = "multiplicative", h = nrow(valid))
ets_model_add <- ets(energy_ts, model = "ZZA")
ets_model_mul <- ets(energy_ts, model = "ZZM")

esm_models <- list(energy_hw_add, energy_hw_mul, ets_model_add, ets_model_mul)

# Summarize all of the ESM models
lapply(esm_models, summary)

# Plot the forecasts of the ESM models
for (model in esm_models) {
    print(
        autoplot(model) +
            autolayer(fitted(model), series = "Fitted") +
            ylab("Megawatts")
    )
}

# Calculate the MAPE and MAE on the validation set for ESM models
calculate_valid_accuracy <- function(model) {
    pred <- forecast::forecast(model, h = nrow(valid))
    valid_error <- valid$mw - pred$mean
    MAPE <- mean(abs(valid_error) / abs(valid$mw)) * 100.0
    MAE <- mean(abs(valid_error))
    model_method <- pred$method
    print(paste("(", model_method, ")", "MAPE:", round(MAPE, 3), "%", "|", "MAE:", round(MAE, 3)))
}

# Loop over forecast metrics for all ESMs
for (model in esm_models) {
    calculate_valid_accuracy(model)
}

# Check for any seasonal unit roots using Canova-Hansen test
# More than 0 so seasonal differencing is needed to address seasonality in our data
energy_ts %>% nsdiffs()

# Take 1 seasonal difference and plot the differenced data
# We notice that there are counter spikes around lag 24--seasonality is multiplicative
# Spikes in the ACF at low lags with tapering PACF indicate non-seasonal MA terms
# Spikes in the PACF at low lags with a tapering ACF indicate possible non-seasonal AR terms
# PACF: Lag 2 shows large spike, quickly drops afterwards. AR(2) potentially
# ACF: Tapering decrease, but we will try AR(2) first without MA terms
# Seasonal: We see spikes as far back as lag 144. Maybe multiple seasonal patterns?
# We will start with simpler models and build up though because the computation cost is high
energy_ts_diff <- energy_ts %>% diff(lag = 24)
energy_ts_diff %>% ggtsdisplay(lag.max = 24 * 6)

# Check for stationarity using ADF test to see if any regular differences needed
# Data seems to be centered around 0 so we use Type I test
# Reject all null hypotheses so we believe data is now stationary
adf.test(energy_ts_diff)

AIC_BIC_list <- list()
training_metrics <- list()

# Creating different sARIMA models
# Checking the residuals
# Collecting their AICs and BICs for model comparison
# Saving the model to a file because computational time grows with complexity
collect_ARIMA_diagnostics <- function(model) {
    model_call <- deparse(summary(model)$call)
    print(model_call)
    resid <- residuals(model)
    resid %>% ggtsdisplay()
    # Need to refer to globally-scoped variables
    AIC_BIC_list[[model_call]] <<- c(AIC = round(summary(model)$aic, 3), BIC = round(summary(model)$bic, 3))
    training_metrics[[model_call]] <<- c(
        MAE = mean(abs(resid)),
        MAPE = mean(abs(resid) / abs(model$x)) * 100.0
    )
    model %>% checkresiduals()
    saveRDS(model, file = paste("models/", model_call, ".rds", sep = ""))
}

# Manual ARIMA model creation

# ARIMA(2, 0, 0)(1, 1, 1)
# AIC 739011.2 | BIC: 739056.9
s_arima_1 <- Arima(energy_ts, order = c(2, 0, 0), seasonal = c(1, 1, 1))
s_arima_1 <- readRDS("models/Arima(y = energy_ts, order = c(2, 0, 0), seasonal = c(1, 1, 1)).rds")
collect_ARIMA_diagnostics(s_arima_1)

# ARIMA(2, 0, 0)(2, 1, 1)
# AIC: 738846.6254 | BIC: 738901.3606
s_arima_2 <- Arima(energy_ts, order = c(2, 0, 0), seasonal = c(2, 1, 1))
s_arima_2 <- readRDS("models/Arima(y = energy_ts, order = c(2, 0, 0), seasonal = c(2, 1, 1)).rds")
collect_ARIMA_diagnostics(s_arima_2)

# ARIMA(2, 0, 0)(2, 1, 2)
# AIC: 738202.382 | BIC: 738266.2408
s_arima_3 <- Arima(energy_ts, order = c(2, 0, 0), seasonal = c(2, 1, 2))
s_arima_3 <- readRDS("models/Arima(y = energy_ts, order = c(2, 0, 0), seasonal = c(2, 1, 2)).rds")
collect_ARIMA_diagnostics(s_arima_3)

# ARIMA(2, 0, 2)(2, 1, 2)
# AIC: 738076.9621 | BIC: 738159.0649
s_arima_4 <- Arima(energy_ts, order = c(2, 0, 2), seasonal = c(2, 1, 2))
s_arima_4 <- readRDS("models/Arima(y = energy_ts, order = c(2, 0, 2), seasonal = c(2, 1, 2)).rds")
collect_ARIMA_diagnostics(s_arima_4)

# ARIMA(2, 0 ,0)(5, 1, 1)
# AIC: 736985.9 | BIC: 737068.0
s_arima_5 <- Arima(energy_ts, order = c(2, 0, 0), seasonal = c(5, 1, 1))
s_arima_5 <- readRDS("models/Arima(y = energy_ts, order = c(2, 0, 0), seasonal = c(5, 1, 1)).rds")
collect_ARIMA_diagnostics(s_arima_5)

# Create an automated sARIMA for comparison

# ARIMA(2, 0, 0)(2, 1, 0) with drift
# AIC: 751522.3 | BIC: 751577.1
auto_s_arima <- auto.arima(energy_ts, seasonal = TRUE)
auto_s_arima <- readRDS("models/auto_sarima.rds")
summary(auto_s_arima)

auto_resid <- residuals(auto_s_arima)
AIC_BIC_list[["Auto sARIMA"]] <- c(AIC = summary(auto_s_arima)$aic, BIC = summary(auto_s_arima)$bic)
training_metrics[["Auto sARIMA"]] <- c(
    MAE = mean(abs(auto_resid)),
    MAPE = mean(abs(auto_resid) / abs(auto_s_arima$x)) * 100.0
)
saveRDS(auto_s_arima, "models/auto_sarima.rds")

# Find top 2 models based on AIC/BIC and goodness-of-fit
# Top 2: ARIMA(2, 0 ,0)(2, 1, 2) and ARIMA(2, 0, 2)(2, 1, 2)
lapply(training_metrics, function(x) sprintf("%.3f", x))
lapply(AIC_BIC_list, function(x) sprintf("%.3f", x))

arima_models <- list(s_arima_3, s_arima_4, s_arima_5)

# Calculate MAPE and MAE for sARIMA forecasts
# MAPE: 5.346 | MAE: 200.479 -- ARIMA(2, 0, 0)(2, 1, 2)
# MAPE: 5.322 | MAE: 199.568 -- ARIMA(2, 0, 2)(2, 1, 2)
# MAPE: 4.792 | MAE: 180.494 -- ARIMA(2, 0, 0)(5, 1, 1)
for (model in arima_models) {
    calculate_valid_accuracy(model)
}

# Plot the actual values and forecasts for the models
# Save the plot as an image in the directory
plot_forecasts <- function(model) {
    pred <- forecast(model, h = nrow(valid))
    df_plot <- data.frame(
        "date" = as.POSIXct(valid$datetime_beginning_ept, format = "%m/%d/%y %H:00"),
        "valid_fc" = valid$mw,
        "pred_fc" = pred$mean
    )
    print(
        ggplot(df_plot, aes(x = date)) +
            geom_line(aes(y = valid_fc, color = "Actual", group = 1), show.legend = TRUE) +
            geom_line(aes(y = pred_fc, color = "Forecast", group = 1), show.legend = TRUE, linetype = "dashed") +
            labs(
                x = "Date",
                y = "Energy Usage (MW)",
                color = "Legend"
            ) +
            scale_color_manual(values = c("Actual" = "blue", "Forecast" = "orange")) +
            scale_x_datetime(date_breaks = "2 days", date_labels = "%m/%d/%y") +
            theme_bw() +
            theme(
                panel.background = element_rect(fill = "#f4f4f4"),
                axis.text.x = element_text(size = 12),
                axis.text.y = element_text(size = 12),
                legend.text = element_text(size = 12),
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(size = 18)
            )
    )
    ggsave(paste(pred$method, ".png", sep = ""), width = 12, height = 6, dpi = 600)
}

for (model in c(arima_models, esm_models)) {
    plot_forecasts(model)
}

valid %>% head()
energy %>% tail()

View(energy %>% filter(str_detect(datetime_beginning_ept, "11/6/16")))

