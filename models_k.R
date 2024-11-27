
''''''''''''''''''''''''''''''''''''''''''''''''''''
#---- SIR model for 2 age groups in 2 settings -----

''''''''''''''''''''''''''''''''''''''''''''''''''''

library(odin)
library(tidyr)
library(ggplot2)

#__________________________________________________________________________

# Note
  # Sw: Susceptible individuals in the community (no chronic disease)
  # Iw: Infected individuals in the community
  # Rw: Recovered individuals in the community
  # Swh, Iwh, Rwh: Corresponding states for individuals in the hospital
  # N: total population
  # beta: Transmission rates
  # gamma: Recovery rates
  # rho: Admissions and discharge rates
#__________________________________________________________________________

# Age group 1
#-------------------------------------------------------

  # Healthy individuals
age1_well_ode <- odin({
  # Derivatives
  deriv(Sw1) <- -beta_w1 * Sw1 * (Iw1/N1) - rho_h1 * Sw1 + rho_wh1 * Swh1
  deriv(Iw1) <- beta_w1 * Sw1 * (Iw1/N1) - gamma_w1 * Iw1 - rho_h1 * Iw1 + rho_wh1 * Iwh1
  deriv(Rw1) <- gamma_w1 * Iw1 - rho_h1 * Rw1 + rho_wh1 * Rwh1

  deriv(Swh1) <- rho_h1 * Sw1 - beta_wh1 * Swh1 * (Iw1 / N1) - rho_wh1 * Swh1
  deriv(Iwh1) <- rho_h1 * Iw1 + beta_wh1 * Swh1 * (Iw1 / N1) - gamma_wh1 * Iwh1 - rho_wh1 * Iwh1
  deriv(Rwh1) <- rho_h1 * Rw1 + gamma_wh1 * Iwh1 - rho_wh1 * Rwh1 
  
  # Initial conditions
  initial(Sw1) <- N1 - Iw1_init - Rw1_init
  initial(Iw1) <- Iw1_init
  initial(Rw1) <- Rw1_init
  initial(Swh1) <- 0
  initial(Iwh1) <- 0
  initial(Rwh1) <- 0
  
  # Parameters and initial values
  Iw1_init <- 10
  Rw1_init <- 0
  N1 <- user(1000)  # Total population for age group 1
  beta_w1 <- user(0.5)  # Transmission rate in the community
  beta_wh1 <- user(1.0)  # Transmission rate in the hospital
  gamma_w1 <- user(0.1)  # Recovery rate in the community
  gamma_wh1 <- user(0.05)  # Recovery rate in the hospital
  rho_h1 <- user(0.02)  # Rate of hospital admission
  rho_wh1 <- user(0.01)  # Rate of discharge from hospital
  
  # Outputs for tracking
  output(S_total) <- Sw1 + Swh1
  output(I_total) <- Iw1 + Iwh1
  output(R_total) <- Rw1 + Rwh1
})

  # Initialize model
  age1_well_mod <- age1_well_ode$new()
  # How long to run
  times <- seq(0,200)
  # Run the model
  pred <- data.frame(age1_well_mod$run(times))

  # Summing compartments to ensure population is conserved
  # total_population <- pred$Sw1 + pred$Iw1 + pred$Rw1 + pred$Swh1 + pred$Iwh1 + pred$Rwh1
  # plot(times, total_population, type = "l", main = "Population Conservation")

  # Reshape data for plotting
  df_plot1 <- pivot_longer(
    pred,
    cols = Sw1:Rwh1,  # Include all relevant compartments
    names_to = "comp",
    values_to = "n"
  )
  
  # Plot using ggplot2
  ggplot(df_plot1, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    scale_color_brewer(palette = "Dark2") +  # Adjust color palette
    labs(
      color = "Compartment",
      y = "Population",
      x = "Time (days)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Check peak infection
  which.max(pred$Iw1) 

  max_infected <- max(pred$Iw1)
  cat("Peak number of infected individuals:", max_infected, "\n")
  cat("Day of peak infection:", which.max(pred$Iw1), "\n")
  
  ggplot(df_plot, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = 17, linetype = "dashed", color = "red") +
    annotate("text", x = 17, y = max_infected, label = "Peak Infection", hjust = -0.1) +
    scale_color_brewer(palette = "Dark2") +
    labs(color = "Compartment", y = "Population", x = "Time (days)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  #-------------------------------------------------------
  # Chronic diseases individuals

  age1_chronic_ode <- odin({
    # Derivatives
    deriv(Sc1) <- -beta_c1 * Sc1 * (Ic1/N1) - rho_c1 * Sc1 + rho_ch1 * Sch1
    deriv(Ic1) <- beta_c1 * Sc1 * (Ic1/N1) - gamma_c1 * Ic1 - rho_c1 * Ic1 + rho_ch1 * Ich1
    deriv(Rc1) <- gamma_c1 * Ic1 - rho_c1 * Rc1 + rho_ch1 * Rch1
    
    deriv(Sch1) <- rho_c1 * Sc1 - beta_ch1 * Sch1 * (Ic1 / N1) - rho_ch1 * Sch1
    deriv(Ich1) <- rho_c1 * Ic1 + beta_ch1 * Sch1 * (Ic1 / N1) - gamma_ch1 * Ich1 - rho_ch1 * Ich1
    deriv(Rch1) <- rho_c1 * Rc1 + gamma_ch1 * Ich1 - rho_ch1 * Rch1 
    
    # Initial conditions
    initial(Sc1) <- N1 - Ic1_init - Rc1_init
    initial(Ic1) <- Ic1_init
    initial(Rc1) <- Rc1_init
    initial(Sch1) <- 0
    initial(Ich1) <- 0
    initial(Rch1) <- 0
    
    # Parameters and initial values
    Ic1_init <- 10
    Rc1_init <- 0
    N1 <- user(1000)  # Total population for age group 1
    beta_c1 <- user(0.5)  # Transmission rate in the community (for chronic patients)
    beta_ch1 <- user(1.2)  # Transmission rate in the hospital (for chronic patients, higher due to vulnerability)
    gamma_c1 <- user(0.05)  # Recovery rate in the community (slower recovery for chronic patients)
    gamma_ch1 <- user(0.03)  # Recovery rate in the hospital (slower recovery in hospital as well)
    rho_c1 <- user(0.05)  # Rate of hospital admission (higher due to the severity of chronic disease)
    rho_ch1 <- user(0.02)  # Rate of discharge from hospital (lower, as recovery is slower)
    
    # Outputs for tracking
    output(S_total) <- Sc1 + Sch1
    output(I_total) <- Ic1 + Ich1
    output(R_total) <- Rc1 + Rch1
  })
  
  # Initialize model
  age1_chronic_mod <- age1_chronic_ode$new()
  
  # Time sequence for the simulation
  times <- seq(0, 200)
  
  # Run the model and store the results in a data frame
  pred <- data.frame(age1_chronic_mod$run(times))
  
  # Reshape data for plotting
  df_plot2 <- pivot_longer(
    pred,
    cols = Sc1:Rch1,  # Include all relevant compartments
    names_to = "comp",
    values_to = "n"
  )
  
  # Plot using ggplot2
  ggplot(df_plot2, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    scale_color_brewer(palette = "Dark2") +  # Adjust color palette
    labs(
      color = "Compartment",
      y = "Population",
      x = "Time (days)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Check peak infection for chronic individuals in the community
  which.max(pred$Ic1)
  
  max_infected <- max(pred$Ic1)
  cat("Peak number of infected chronic individuals:", max_infected, "\n")
  cat("Day of peak infection:", which.max(pred$Ic1), "\n")
  
  # Optional: Adding a peak infection marker on the plot
  ggplot(df_plot1, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    geom_vline(xintercept = 30, linetype = "dashed", color = "red") +
    annotate("text", x = 30, y = max_infected, label = "Peak Infection", hjust = -0.1) +
    scale_color_brewer(palette = "Dark2") +
    labs(color = "Compartment", y = "Population", x = "Time (days)") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
#__________________________________________________________________________
  
# Age group 2
#------------------------------------------------------- 
  
  # Healthy individuals
  age2_well_ode <- odin({
    # Derivatives
    deriv(Vw2) <- vw2 * Sw2
    deriv(Sw2) <- -vw2 * Sw2 -beta_w2 * Sw2 * (Iw2/N2) - rho_w2 * Sw2 + rho_wh2 * Swh2
    deriv(Iw2) <- beta_w2 * Sw2 * (Iw2 / N2) - gamma_w2 * Iw2 - rho_w2 * Iw2 + rho_wh2 * Iwh2
    deriv(Rw2) <- gamma_w2 * Iw2 - rho_w2 * Rw2 + rho_wh2 * Rwh2
    
    deriv(Vwh2) <- vwh2 * Swh2
    deriv(Swh2) <- -vwh2 * Swh2 + rho_w2 * Sw2 - beta_wh2 * Swh2 * (Iw2 / N2) - rho_wh2 * Swh2
    deriv(Iwh2) <- rho_w2 * Iw2 + beta_wh2 * Swh2 * (Iw2 / N2) - gamma_wh2 * Iwh2 - rho_wh2 * Iwh2
    deriv(Rwh2) <- rho_w2 * Rw2 + gamma_wh2 * Iwh2 - rho_wh2 * Rwh2
    
    # Initial conditions for Age Group 2
    initial(Sw2) <- N2 - Iw2_init - Rw2_init
    initial(Iw2) <- Iw2_init
    initial(Rw2) <- Rw2_init
    initial(Swh2) <- 0
    initial(Iwh2) <- 0
    initial(Rwh2) <- 0
    initial(Vw2) <- 0
    initial(Vwh2) <- 0
    
    # Parameters for Age Group 2
    Iw2_init <- 10         # Initial number of infected individuals
    Rw2_init <- 0           # Initial number of recovered individuals
    N2 <- user(1000)       # Total population for age group 2
    vw2 <- user(0.3)       # Vaccination rate of healthy individuals in the community (e.g., 30% of the susceptible population gets vaccinated per time step)
    vwh2 <- user(0.5)      # Vaccination rate of healthy individuals in the hospital
    beta_w2 <- user(0.3)   # Transmission rate in the community (susceptible -> infected)
    beta_wh2 <- user(0.6)  # Transmission rate in the hospital (susceptible -> infected)
    gamma_w2 <- user(0.2)  # Recovery rate in the community
    gamma_wh2 <- user(0.1)# Recovery rate in the hospital
    rho_w2 <- user(0.02)   # Rate of hospital admission
    rho_wh2 <- user(0.01)  # Rate of discharge from hospital
    
    # Outputs for tracking
    output(S_total) <- Sw2 + Swh2
    output(I_total) <- Iw2 + Iwh2
    output(R_total) <- Rw2 + Rwh2
    output(V_total) <- Vw2 + Vwh2
  })
    
  # Initialize model
  age2_well_mod <- age2_well_ode$new()
  # How long to run
  times <- seq(0,200)
  # Run the model
  pred <- data.frame(age2_well_mod$run(times))
  
  # Summing compartments to ensure population is conserved
  # total_population <- pred$Sw1 + pred$Iw1 + pred$Rw1 + pred$Swh1 + pred$Iwh1 + pred$Rwh1
  # plot(times, total_population, type = "l", main = "Population Conservation")
  
  # Reshape data for plotting
  df_plot3 <- pivot_longer(
    pred,
    cols = Sw2:Rwh2,  # Include all relevant compartments
    names_to = "comp",
    values_to = "n"
  )
  
  # Plot using ggplot2
  ggplot(df_plot3, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    scale_color_brewer(palette = "Dark2") +  # Adjust color palette
    labs(
      color = "Compartment",
      y = "Population",
      x = "Time (days)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Check peak infection
  which.max(pred$Iw2) 

#------------------------------------------------------- 
  
  # Chronic diseases individuals 
  age2_chronic_ode <- odin({
    # Derivatives for Age Group 2 (Chronic Patients)
    deriv(Vc2) <- vc2 * Sc2
    deriv(Sc2) <- -vc2 * Sc2 - beta_c2 * Sc2 * (Ic2 / N2) - rho_c2 * Sc2 + rho_ch2 * Sch2
    deriv(Ic2) <- beta_c2 * Sc2 * (Ic2 / N2) - gamma_c2 * Ic2 - rho_c2 * Ic2 + rho_ch2 * Ich2
    deriv(Rc2) <- gamma_c2 * Ic2 - rho_c2 * Rc2 + rho_ch2 * Rch2
    
    deriv(Vch2) <- vch2 * Sch2
    deriv(Sch2) <- -vch2 * Sch2 + rho_c2 * Sc2 - beta_ch2 * Sch2 * (Ic2 / N2) - rho_ch2 * Sch2
    deriv(Ich2) <- rho_c2 * Ic2 + beta_ch2 * Sch2 * (Ic2 / N2) - gamma_ch2 * Ich2 - rho_ch2 * Ich2
    deriv(Rch2) <- rho_c2 * Rc2 + gamma_ch2 * Ich2 - rho_ch2 * Rch2
    
    # Initial conditions for Age Group 2 (Chronic Patients)
    initial(Sc2) <- N2 - Ic2_init - Rc2_init
    initial(Ic2) <- Ic2_init
    initial(Rc2) <- Rc2_init
    initial(Sch2) <- 0
    initial(Ich2) <- 0
    initial(Rch2) <- 0
    initial(Vc2) <- 0
    initial(Vch2) <- 0
    
    # Parameters for Age Group 2 (Chronic Patients)
    Ic2_init <- 10        # Initial number of infected chronic individuals
    Rc2_init <- 0         # Initial number of recovered chronic individuals
    N2 <- user(1000)      # Total population for age group 2
    vc2 <- user(0.1)      # Vaccination rate of chronic individuals in the community (e.g., 10% of the susceptible population gets vaccinated per time step)
    vch2 <- user(0.15)    # Vaccination rate of chronic individuals in the hospital (e.g., 15% of the susceptible population gets vaccinated per time step)
    beta_c2 <- user(0.2)  # Transmission rate in the community (susceptible -> infected) for chronic individuals
    beta_ch2 <- user(0.4) # Transmission rate in the hospital (susceptible -> infected) for chronic individuals
    gamma_c2 <- user(0.1) # Recovery rate in the community for chronic individuals
    gamma_ch2 <- user(0.05) # Recovery rate in the hospital for chronic individuals
    rho_c2 <- user(0.05)  # Rate of hospital admission for chronic individuals (higher than healthy)
    rho_ch2 <- user(0.02) # Rate of discharge from hospital for chronic individuals (lower than healthy)
    
    # Outputs for tracking
    output(S_total) <- Sc2 + Sch2
    output(I_total) <- Ic2 + Ich2
    output(R_total) <- Rc2 + Rch2
    output(V_total) <- Vc2 + Vch2
  })
  
  # Initialize model
  age2_chronic_mod <- age2_chronic_ode$new()
  # How long to run
  times <- seq(0,200)
  # Run the model
  pred <- data.frame(age2_chronic_mod$run(times))
  
  # Reshape data for plotting
  df_plot4 <- pivot_longer(
    pred,
    cols = Sc2:Rch2,  # Include all relevant compartments
    names_to = "comp",
    values_to = "n"
  )
  
  # Plot using ggplot2
  ggplot(df_plot4, aes(x = t, y = n, color = comp)) +
    geom_line(linewidth = 1.2) +
    scale_color_brewer(palette = "Dark2") +  # Adjust color palette
    labs(
      color = "Compartment",
      y = "Population",
      x = "Time (days)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Check peak infection
  which.max(pred$Ic2) 
  













