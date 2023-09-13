#////////////
# Simulates times series, fits basic OLS linear model to them, and then does AIC model selection
#///////////

library(MASS)

# make a list to hold all of the simulated datasets and their fitted models


# Step 1: Simulate time series --------------------------------------------
# use basic_timeseries() function (no seasonality, 1 sample/time step)
#simDat <- basic_timeseries(K = 50, num_strong = 5, n = 60, freq = 1, prob_cycle = 0)

ts_AICselection <- lapply(as.list(1:1000),
                          FUN = function(x)
                            basic_timeseries(K = 50, num_strong = 5, n = 60, freq = 1, prob_cycle = 0)
)


for (i in 1:length(ts_AICselection)) {
  # put data into a data.frame
  simDat_df <- data.frame("y" = ts_AICselection[[i]]$y, ts_AICselection[[i]]$X[,2:51])
  # make names better
  names(simDat_df) <- c("y",paste0("driver_",2:51))
  # fit full global model
  mod <- lm(  y ~   driver_2 +   driver_3 +   driver_4+   driver_5+   driver_6+   driver_7+   driver_8+   driver_9
              +   driver_10  +   driver_11   +   driver_12   +   driver_13  +   driver_14  +   driver_15  +   driver_16  +   driver_17  +   driver_18  +   driver_19
              +   driver_20  +   driver_21   +   driver_22   +   driver_23  +   driver_24  +   driver_25  +   driver_26  +   driver_27  +   driver_28  +   driver_29
              +   driver_30  +   driver_31   +   driver_32   +   driver_33  +   driver_34  +   driver_35  +   driver_36  +   driver_37  +   driver_38  +   driver_39
              +   driver_40  +   driver_41   +   driver_42   +   driver_43  +   driver_44  +   driver_45  +   driver_46  +   driver_47  +   driver_48  +   driver_49
              +   driver_50 +   driver_51, data = simDat_df)

  # fit null model
  mod_null <- lm(y ~ 1, data = simDat_df)
  # Do AIC model selection
  stepAIC <- MASS::stepAIC(mod_null, direction = "forward", scope = list(lower = mod_null,
                                                                         upper = mod))

  ## save the output
  # save the full model
  ts_AICselection[[i]]$globalModel <- mod
  # save the AIC selected model
  ts_AICselection[[i]]$stepAIC_results<- stepAIC
}
