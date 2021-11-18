library(shiny)
library(shinyWidgets)
library(tidyverse)
library(mvtnorm)
library(patchwork)

ui <- fluidPage(withMathJax(), #Allows Latex text to be used
  tags$div(HTML("<script type='text/x-mathjax-config'>
  MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
  });  </script> ")),
  
  titlePanel("Using the Lkj Prior"),
  
  "In this demonstration, we generate data from a simple hiererchical model with both a random intercept and random slope.", 
  "We then use the brms package to fit a Bayesian model with an Lkj prior.", 
  "The goal is to demonstrate how the choice of parameter for the Lkj prior effects our posterior inference.",
  
  h3("Model Form"),
  
  "We begin by describing the model's form.",
  "Let $j \\in \\{1, \\ldots, J\\}$ index the groups and $i \\in \\{1, \\ldots, n_j\\}$ index observations within each group.",
  "Let $y_{ij}$ be the response variable and $x_{ij}$ be a predictor.", 
  "Then we generate data and model via the mixed effects model below.",
  
  "$$ y_{ij} = (\\beta_{0} + b_{0j}) + (\\beta_{1} + b_{1j}) x_{ij} + \\varepsilon_{ij}, \\qquad \\varepsilon_{ij} \\overset{iid}{\\sim} \\mathcal{N}(0, \\sigma^2), \\qquad",
  "\\begin{pmatrix} b_{0j} \\\\ b_{1j} \\end{pmatrix} \\overset{iid}{\\sim} \\mathcal{N}_2(\\mathbf{0}, \\mathbf{\\Sigma}) $$",
  
  "Decompose $\\mathbf{\\Sigma}$ into", 
  "$$ \\mathbf{\\Sigma} = \\begin{pmatrix} \\tau_0 & 0 \\\\ 0 & \\tau_1 \\end{pmatrix}\\mathbf{\\Omega}\\begin{pmatrix} \\tau_0 & 0 \\\\ 0 & \\tau_1 \\end{pmatrix},$$",
  "for correlation matrix $\\mathbf{\\Omega} = \\begin{pmatrix} 1 & \\rho \\\\ \\rho & 1 \\end{pmatrix}$.",
  "The amount of correlation between the random intercept and slope can be selected below.",
  
  h3("Generate Data"),
  
  sliderTextInput(inputId = "rho", label = withMathJax(HTML("Choose A True Value for $\\rho$")),
                  selected = 0, choices = c(seq(-0.95, -0.15, 0.05), seq(-0.1, 0.95, 0.05))),
  
  "The other counts and parameters are set to $J = 8$, $n_j = 25$ for all $j$, $\\beta_0 = 5$, $\\beta_1 = 2$, $\\tau_0 = 2$, $\\tau_1 = 1$, and $\\sigma = 1$.",
  
  "Data generated with this $\\rho$ is plotted here.",
  plotOutput("data_plot"),
  
  h3("Prior Structure"),
  
  "For simplicity, we place weakly informative priors on $\\beta_0$, $\\beta_1$, $\\tau_0$, and $\\tau_1$. The prior for $\\mathbf{\\Omega}$ is",
  "$$ \\mathbf{\\Omega} \\sim LkjCorr(\\eta) $$",
  "Choose the prior parameter $\\eta$ below and note the induced prior on the correlation, $\\rho$.",
  sliderTextInput(inputId = "eta", label = withMathJax(HTML("Choose The Prior Parameter $\\eta$")),
                  selected = 1, choices = c(0.1, seq(0.2, 0.8, 0.2), 1:5, 10)),
  plotOutput("prior_plot"),
  
  h3("Posterior"),
  
  "The posterior is computed using the brms package with the priors above.",
  "The prior is the red density, the posterior is the blue density, and the true $\\rho$ is the dashed line.",
  "The sliders to select values for $\\rho$ and $\\eta$ are reporoduced below to make it easy to try various combinations of the two.",
  
  sliderTextInput(inputId = "rho2", label = withMathJax(HTML("Choose A True Value for $\\rho$")),
                  selected = 0, choices = c(seq(-0.95, -0.15, 0.05), seq(-0.1, 0.95, 0.05))),
  sliderTextInput(inputId = "eta2", label = withMathJax(HTML("Choose The Prior Parameter $\\eta$")),
                  selected = 1, choices = c(0.1, seq(0.2, 0.8, 0.2), 1:5, 10)),
  
  plotOutput("posterior_plot")
  )

server <- function(input, output, session){
  #Make the two rho and eta sliders match
  observeEvent(input$eta2, {
    updateSliderTextInput(session = session, inputId = "eta", 
                          selected = input$eta2)
  })
  observeEvent(input$eta, {
    updateSliderTextInput(session = session, inputId = "eta2", 
                          selected = input$eta)
  })
  observeEvent(input$rho2, {
    updateSliderTextInput(session = session, inputId = "rho", 
                          selected = input$rho2)
  })
  observeEvent(input$rho, {
    updateSliderTextInput(session = session, inputId = "rho2", 
                          selected = input$rho)
  })
  
  #Generate the data (with a pre-set seed consistent with cached results)
  data <- reactive({
    J <- 8
    n <- rep(25,J)
    beta0 <- 5
    beta1 <- 2
    tau0 <- 2
    tau1 <- 1
    sigma <- 1
    
    tau_mat <- diag(c(tau0, tau1))
    Omega <- matrix(c(1,input$rho,input$rho,1),ncol=2)
    Sigma <- tau_mat %*% Omega %*% tau_mat
    
    set.seed(1)
    beta <- rmvnorm(J, mean = c(beta0, beta1),  sigma = Sigma) #Sample coefficients
    
    x <- c(); y <- c()
    for(j in 1:J){
      xj <- rnorm(n[j])
      yj <- rnorm(n[j], beta[j,1] + beta[j,2]*xj,sigma)
      x <- c(x,xj); y <- c(y,yj)
    }
    
    data.frame(x, y, group = rep(1:J, n), 
               Intercept = rep(beta[,1], n),
               Slope = rep(beta[,2], n))})
  
  #Plot the data and the generated ("true") intercepts/slopes
  output$data_plot <- renderPlot({
    p1 <- ggplot(data(), aes(x,y, color = factor(group))) + geom_point() + 
      labs(x = expression("x"[ij]), y = expression("y"[ij]), color = "Group j",
           title = "Observed Data") +
      theme_classic()
    
    p2 <- ggplot(data(), aes(x = Intercept, y = Slope, color = factor(group))) +
      geom_point() + theme_classic() +
      labs(x = expression(beta[0]*" + b"[0*"j"]), 
           y = expression(beta[1]*" + b"[1*"j"]), 
           color = "Group j", 
           title = "True Intercepts and Slopes for Each Group")
    
    p1 + p2
    })
  
  #Plot the prior for a given choice of eta
  output$prior_plot <- renderPlot({
    rho_dens <- function(rho){0.5*dbeta((rho+1)/2, input$eta, input$eta)}
    
    ggplot(data.frame(x=c(-1,1)), aes(x)) + 
      geom_function(fun = rho_dens, color = "red") +
      stat_function(fun = rho_dens, geom = "area", fill = "red", alpha = 0.1,
                    aes(color = "Prior")) +
      geom_vline(xintercept = input$rho, linetype = "dashed") + 
      labs(x = expression(rho), color = "") +
      theme_classic()
  })
  
  #Read in the cached posterior samples for given rho/eta combination
  post_samp <- reactive({
    ############# CHANGE THIS FILE PATH, IF NECCESARY #######################
    read_csv(paste0("hw4_cache/",
    #########################################################################
                    "eta_", as.character(input$eta), 
                    "_rho_", as.character(input$rho), ".csv"),
             col_types = "d")
  })
  
  #Plot the prior and posterior for a given rho/eta combination
  output$posterior_plot <- renderPlot({
    rho_dens <- function(rho){0.5*dbeta((rho+1)/2, input$eta, input$eta)}

     ggplot(data.frame(post = pull(post_samp())), aes(x = post)) +
      geom_function(fun = rho_dens, color = "red") +
       stat_function(fun = rho_dens, geom = "area", fill = "red", alpha = 0.1,
                     aes(color = "Prior")) +
       geom_density(alpha = 0.1, color = "blue", fill = "blue",
                    aes(linetype = "Posterior")) +
       geom_vline(xintercept = input$rho, linetype = "dashed") +
       xlim(-1,1) +
       labs(x = expression(rho), linetype = "", color = "") +
       theme_classic()
  })
}

shinyApp(ui = ui, server = server)
