library(tidyverse, warn.conflicts = FALSE)
library(forecast)
library(tseries)

# Read in train and valid data
train <- read_csv("data/clean_train.csv")
valid <- read_csv("data/clean_valid.csv")
test <- read_csv("data/clean_valid_2.csv")

season_length <- 24
train_ts <- ts(train$mw, frequency = season_length)
valid_ts <- ts(valid$mw, frequency = season_length)
test_ts <- ts(test$mw, frequency = season_length)

# Test for seasonal differences for neural net model later
train_ts %>% nsdiffs()

# Take a difference of 24 lags to make data stationary
train_diff <- train_ts %>% diff(season_length)

# Plot differenced data
autoplot(train_diff)

# Checking for any regular differences to be taken
train_diff %>% ndiffs()

set.seed(1234)

# Display the correlation plots of differenced data to get an idea of which AR terms to use
train_ts %>%
  diff(lag = season_length) %>%
  ggtsdisplay()

# Fit a neural net model with p = 3 and P = 3
nn_model <- nnetar(diff(train_ts, season_length), size = 10, p = 3, P = 3)

# Forecast the next 168 hours using the neural net model
nn_forecast <- forecast::forecast(nn_model, h = nrow(valid))

# Adjust the forecasted values to be on the same scale as the original data based off the last season in the training data
pass_forecast <- train_ts[(length(train_ts) - season_length + 1):length(train_ts)] + nn_forecast$mean[1:24]

# Adjust the forecasted values to be on the same scale as the original data
for (i in 25:168) {
  pass_forecast[i] <- pass_forecast[i - season_length] + nn_forecast$mean[i]
}
pass_forecast <- ts(pass_forecast, frequency = season_length)

# Calculate the MAPE and MAE of the forecasted values
# [Valid1 Metrics] MAE: 335.875630 | MAPE: 9.021000 
nn_error <- valid_ts - pass_forecast
nn_MAE <- mean(abs(nn_error))
nn_MAPE <- mean(abs(nn_error) / abs(valid_ts)) * 100.0
sprintf("MAE: %f", nn_MAE)
sprintf("MAPE: %f", nn_MAPE)

# Plot the forecasted values vs. actual values in validation data
ggplot(valid, aes(datetime_beginning_ept)) +
  geom_line(aes(y = mw, color = "Actual"), show.legend = TRUE) +
  geom_line(aes(y = pass_forecast, color = "Forecast"), show.legend = TRUE, linetype = "dashed") +
  labs(
    x = "Date",
    y = "Energy Usage (MW)"
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
ggsave("neural-net-forecast.png", width = 12, height = 6, dpi = 600)

# Combine train and validation into one dataset
combined_df <- bind_rows(train, valid)
combined_ts <- ts(combined_df$mw, frequency = season_length)

# Fit the final model with same parameters as before
nn_final_model <- nnetar(diff(combined_ts, season_length), size = 10, p = 3, P = 3)

# Forecast the next 168 hours using the final model
nn_final_forecast <- forecast::forecast(nn_final_model, h = nrow(test))

# Adjust the forecasted values to be on the same scale as the original data based off the last season in the combined data
pass_test_forecast <- combined_ts[(length(combined_ts) - season_length + 1):length(combined_ts)] + nn_final_forecast$mean[1:24]

# Adjust the forecasted values to be on the same scale as the original data
for (i in 25:nrow(test)) {
  pass_test_forecast[i] <- pass_test_forecast[i - season_length] + nn_final_forecast$mean[i]
}

# Calculate the MAPE and MAE of the forecasted values
# [Test Metrics] MAE: 214.194150 | MAPE: 5.678016
nn_test_error <- test_ts - pass_test_forecast
nn_test_MAPE <- mean(abs(nn_test_error) / abs(test$mw)) * 100.0
nn_test_MAE <- mean(abs(nn_test_error))
sprintf("Test MAE: %f", nn_test_MAE)
sprintf("Test MAPE: %f", nn_test_MAPE)