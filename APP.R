library(bs4Dash)
library(survival)
library(survminer)
library(ggplot2)
library(randomForestSRC)
library(glmnet)
library(rms)
library(DT)
library(reshape2)

# Helper: Create differentiated newdata for predictions.
# For numeric covariates, subtract 5 for the first treatment level and add 5 for the second.
getDifferentiatedData <- function(df, covariates) {
  if (!("treatment" %in% names(df))) {
    return(NULL)
  }
  trtLevels <- levels(df$treatment)
  newdata_list <- lapply(trtLevels, function(trt) {
    newrow <- data.frame(treatment = trt, stringsAsFactors = FALSE)
    for (v in covariates) {
      if (is.numeric(df[[v]])) {
        med <- median(df[[v]], na.rm = TRUE)
        delta <- ifelse(trt == trtLevels[1], -5, 5)
        newrow[[v]] <- med + delta
      } else {
        levs <- levels(df[[v]])
        newrow[[v]] <- if (length(levs) >= 2) {
          if (trt == trtLevels[1]) levs[1] else levs[2]
        } else { 
          levs[1] 
        }
      }
    }
    newrow
  })
  newdata <- do.call(rbind, newdata_list)
  newdata$treatment <- factor(newdata$treatment, levels = trtLevels)
  newdata
}

# Create sample data with an artificial treatment effect.
set.seed(123)
n <- 500
sampleData <- data.frame(
  time = rexp(n, rate = 0.05),
  status = rbinom(n, 1, prob = 0.7),
  treatment = factor(sample(c("CABG", "PCI"), n, replace = TRUE)),
  age = round(runif(n, 50, 80)),
  diabetes = factor(sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.3, 0.7))),
  gender = factor(sample(c("Male", "Female"), n, replace = TRUE)),
  hypertension = factor(sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.6, 0.4))),
  smoking = factor(sample(c("Yes", "No"), n, replace = TRUE, prob = c(0.4, 0.6))),
  num_vessels = sample(1:3, n, replace = TRUE)
)
# Adjust time for treatment effect.
sampleData$time <- ifelse(sampleData$treatment == "CABG",
                          sampleData$time * 0.8,
                          sampleData$time * 1.2)

# UI using bs4Dash
ui <- dashboardPage(
  title = "Advanced IPD Meta-Analysis Dashboard",
  header = dashboardHeader(title = "Advanced IPD Meta-Analysis"),
  sidebar = dashboardSidebar(
    sidebarMenu(
      menuItem("Data", tabName = "data", icon = icon("table")),
      menuItem("Analysis", tabName = "analysis", icon = icon("chart-line")),
      menuItem("Diagnostics", tabName = "diagnostics", icon = icon("stethoscope")),
      menuItem("Risk Prediction", tabName = "risk", icon = icon("user-md")),
      menuItem("Simulation", tabName = "simulation", icon = icon("flask")),
      menuItem("Validation", tabName = "validation", icon = icon("check-circle")),
      menuItem("Extra", tabName = "extra", icon = icon("cogs")),
      menuItem("Advanced", tabName = "advanced", icon = icon("wrench"))
    )
  ),
  body = dashboardBody(
    bs4TabItems(
      # (1) Data Tab
      bs4TabItem(tabName = "data",
                 fluidRow(
                   box(width = 12,
                       title = "Upload or Use Sample Data",
                       fileInput("datafile", "Upload CSV File",
                                 accept = c("text/csv", "text/plain", ".csv")),
                       downloadButton("downloadData", "Download Data (CSV)"),
                       DTOutput("dataTable")
                   )
                 )
      ),
      # (2) Analysis Tab
      bs4TabItem(tabName = "analysis",
                 fluidRow(
                   box(width = 4,
                       title = "Model Options",
                       selectInput("modelType", "Choose Model Type",
                                   choices = c("Cox Proportional Hazards",
                                               "Cox with Splines",
                                               "Time-Dependent Cox",
                                               "Random Survival Forest",
                                               "Penalized Cox (Lasso)")),
                       checkboxGroupInput("covariates", "Select Covariates",
                                          choices = c("age", "diabetes", "gender", "hypertension", "smoking", "num_vessels"),
                                          selected = c("age", "diabetes")),
                       numericInput("alpha", "Significance Level (α)", value = 0.05, min = 0.01, max = 0.1, step = 0.01),
                       helpText("Requires a column named 'treatment' with levels 'CABG' and 'PCI'.")
                   ),
                   box(width = 8,
                       title = "Model Summary",
                       verbatimTextOutput("modelSummary")
                   )
                 ),
                 fluidRow(
                   box(width = 12,
                       title = "Survival Plot",
                       plotOutput("survPlot")
                   )
                 )
      ),
      # (3) Diagnostics Tab
      bs4TabItem(tabName = "diagnostics",
                 fluidRow(
                   box(width = 12,
                       title = "Model Diagnostics / Variable Importance",
                       plotOutput("diagnosticPlot")
                   )
                 )
      ),
      # (4) Risk Prediction Tab
      bs4TabItem(tabName = "risk",
                 fluidRow(
                   box(width = 12,
                       title = "Risk Prediction (Cox PH Model)",
                       plotOutput("riskPlot")
                   )
                 )
      ),
      # (5) Simulation Tab
      bs4TabItem(tabName = "simulation",
                 fluidRow(
                   box(width = 4,
                       title = "Simulation Settings",
                       sliderInput("simAge", "Simulated Age", min = 40, max = 90, value = 65),
                       radioButtons("simTreatment", "Treatment", choices = c("CABG", "PCI"), selected = "CABG")
                   ),
                   box(width = 8,
                       title = "Simulated Survival Curves (Cox PH Model)",
                       plotOutput("simulationPlot")
                   )
                 )
      ),
      # (6) Validation Tab
      bs4TabItem(tabName = "validation",
                 fluidRow(
                   box(width = 4,
                       title = "Validation Settings",
                       numericInput("nBootstrap", "Number of Bootstrap Replicates", value = 100, min = 10, max = 500),
                       actionButton("runBootstrap", "Run Bootstrap Validation")
                   ),
                   box(width = 8,
                       title = "Bootstrap Validation Results (C-index)",
                       verbatimTextOutput("validationResults")
                   )
                 )
      ),
      # (7) Extra Tab: Additional Plots and Text Output Options
      bs4TabItem(tabName = "extra",
                 fluidRow(
                   box(width = 12,
                       title = "Extra Analysis Options",
                       radioButtons("extraOutputType", "Select Extra Output:",
                                    choices = c("Age Histogram", "Time Boxplot", "Dataset Summary"),
                                    selected = "Age Histogram", inline = TRUE),
                       uiOutput("extraUI")
                   )
                 )
      ),
      # (8) Advanced Analysis Tab: Scatterplot Matrix & Correlation Heatmap
      bs4TabItem(tabName = "advanced",
                 fluidRow(
                   box(width = 6,
                       title = "Scatterplot Matrix",
                       plotOutput("scatterPlotMatrix")
                   ),
                   box(width = 6,
                       title = "Correlation Heatmap",
                       plotOutput("corHeatmap")
                   )
                 )
      )
    )
  )
)

# SERVER
server <- function(input, output, session) {
  
  ## Reactive Data
  reactiveData <- reactive({
    df <- if (is.null(input$datafile)) {
      sampleData
    } else {
      read.csv(input$datafile$datapath, stringsAsFactors = TRUE)
    }
    if ("treatment" %in% names(df)) {
      df$treatment <- as.factor(df$treatment)
    }
    df
  })
  
  output$dataTable <- renderDT({
    DT::datatable(reactiveData())
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { "ipd_data.csv" },
    content = function(file) {
      write.csv(reactiveData(), file, row.names = FALSE)
    }
  )
  
  ## Build Model Formula
  modelFormula <- reactive({
    covs <- input$covariates
    formText <- paste("Surv(time, status) ~ treatment",
                      if (length(covs) > 0) paste("+", paste(covs, collapse = " + ")) else "")
    as.formula(formText)
  })
  
  ## Fit Model
  fitModel <- reactive({
    df <- reactiveData()
    modType <- input$modelType
    if (modType == "Cox Proportional Hazards") {
      fit <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
    } else if (modType == "Cox with Splines") {
      if ("age" %in% input$covariates) {
        other_covs <- setdiff(input$covariates, "age")
        formText <- paste("Surv(time, status) ~ treatment + pspline(age)",
                          if (length(other_covs) > 0) paste("+", paste(other_covs, collapse = " + ")) else "")
        fit <- coxph(as.formula(formText), data = df, x = TRUE, y = TRUE)
      } else {
        fit <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
      }
    } else if (modType == "Time-Dependent Cox") {
      if ("age" %in% input$covariates) {
        other_covs <- setdiff(input$covariates, "age")
        formText <- paste("Surv(time, status) ~ treatment + tt(age)",
                          if (length(other_covs) > 0) paste("+", paste(other_covs, collapse = " + ")) else "")
        fit <- coxph(as.formula(formText), data = df, x = TRUE, y = TRUE,
                     tt = function(x, t, ...) { x * log(t + 1) })
      } else {
        fit <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
      }
    } else if (modType == "Random Survival Forest") {
      covNames <- c("treatment", input$covariates)
      fit <- rfsrc(Surv(time, status) ~ ., data = df[, c("time", "status", covNames)])
    } else if (modType == "Penalized Cox (Lasso)") {
      x <- model.matrix(modelFormula(), data = df)[, -1]
      y <- with(df, Surv(time, status))
      fit <- cv.glmnet(x, y, family = "cox")
    }
    fit
  })
  
  output$modelSummary <- renderPrint({
    fit <- fitModel()
    modType <- input$modelType
    if (modType %in% c("Cox Proportional Hazards", "Cox with Splines", "Time-Dependent Cox")) {
      summary(fit)
    } else if (modType == "Random Survival Forest") {
      print(fit)
    } else if (modType == "Penalized Cox (Lasso)") {
      cat("Optimal Lambda:", fit$lambda.min, "\n")
      print(coef(fit, s = "lambda.min"))
    }
  })
  
  ## Kaplan-Meier Survival Plot
  output$survPlot <- renderPlot({
    df <- reactiveData()
    modType <- input$modelType
    fit <- fitModel()
    if (modType == "Time-Dependent Cox") {
      fit <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
    }
    if (modType %in% c("Cox Proportional Hazards", "Cox with Splines", "Time-Dependent Cox")) {
      newdata <- getDifferentiatedData(df, input$covariates)
      surv_fit <- tryCatch(survfit(fit, newdata = newdata),
                           error = function(e) { return(NULL) })
      if (is.null(surv_fit) || (is.null(surv_fit$strata) && length(surv_fit$time) == 0)) {
        plot.new()
        text(0.5, 0.5, "No survival curves available (null model).")
      } else {
        legendLabels <- if (nrow(newdata) == 2) {
          levels(df$treatment)
        } else {
          "Overall"
        }
        ggsurv <- tryCatch({
          ggsurvplot(surv_fit, data = df, risk.table = TRUE, pval = TRUE,
                     legend.title = "Treatment", legend.labs = legendLabels)
        }, error = function(e) { return(NULL) })
        if (is.null(ggsurv)) {
          plot.new()
          text(0.5, 0.5, "ggsurvplot failed.")
        } else {
          print(ggsurv$plot)
        }
      }
    } else if (modType == "Random Survival Forest") {
      plot(fit$chf, xlab = "Time", ylab = "Cumulative Hazard",
           main = "RSF Cumulative Hazard Functions")
    } else if (modType == "Penalized Cox (Lasso)") {
      x <- model.matrix(modelFormula(), data = df)[, -1]
      risk <- predict(fit, newx = x, s = "lambda.min", type = "link")
      df$riskGroup <- "Overall"
      km_fit <- survfit(Surv(time, status) ~ riskGroup, data = df)
      ggsurv <- tryCatch({
        ggsurvplot(km_fit, data = df, risk.table = TRUE, pval = TRUE,
                   legend.title = "Risk Group", legend.labs = "Overall")
      }, error = function(e) { return(NULL) })
      if (is.null(ggsurv)) {
        plot.new()
        text(0.5, 0.5, "Penalized Cox survival plot failed.")
      } else {
        print(ggsurv$plot)
      }
    }
  })
  
  ## Diagnostic Plot Output
  output$diagnosticPlot <- renderPlot({
    df <- reactiveData()
    modType <- input$modelType
    fit <- fitModel()
    if (modType %in% c("Cox Proportional Hazards", "Cox with Splines")) {
      czph <- cox.zph(fit)
      plot(czph, main = "Schoenfeld Residuals")
    } else if (modType == "Time-Dependent Cox") {
      fit2 <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
      czph <- cox.zph(fit2)
      plot(czph, main = "Schoenfeld Residuals (Cox PH)")
    } else if (modType == "Random Survival Forest") {
      vi <- fit$importance
      if (!is.null(vi)) {
        barplot(vi, main = "RSF Variable Importance", col = "skyblue", las = 2)
      } else {
        plot.new()
        text(0.5, 0.5, "No variable importance available for RSF.")
      }
    } else if (modType == "Penalized Cox (Lasso)") {
      plot.new()
      text(0.5, 0.5, "Diagnostics not available for Penalized Cox (Lasso).")
    }
  })
  
  ## Risk Prediction Plot (Cox PH Model)
  output$riskPlot <- renderPlot({
    df <- reactiveData()
    fit_cox <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
    newdata <- getDifferentiatedData(df, input$covariates)
    newdata_pred <- newdata[1, , drop = FALSE]
    surv_fit <- tryCatch(survfit(fit_cox, newdata = newdata_pred),
                         error = function(e) { return(NULL) })
    if (is.null(surv_fit)) {
      plot.new()
      text(0.5, 0.5, "Risk prediction not available (null model).")
    } else {
      lab <- unique(newdata_pred$treatment)
      ggsurv <- tryCatch({
        ggsurvplot(surv_fit, data = df, risk.table = TRUE, pval = FALSE,
                   legend.title = "Treatment", legend.labs = lab,
                   title = "Predicted Survival (Cox PH Model)")
      }, error = function(e) { return(NULL) })
      if (is.null(ggsurv)) {
        plot.new()
        text(0.5, 0.5, "Risk prediction plot failed.")
      } else {
        print(ggsurv$plot)
      }
    }
  })
  
  ## Advanced Tab: Scatterplot Matrix & Correlation Heatmap
  output$scatterPlotMatrix <- renderPlot({
    df <- reactiveData()
    nums <- df[sapply(df, is.numeric)]
    if (ncol(nums) > 1) {
      pairs(nums, main = "Scatterplot Matrix of Numeric Variables", col = "blue")
    } else {
      plot.new()
      text(0.5, 0.5, "Not enough numeric variables for scatterplot matrix.")
    }
  })
  
  output$corHeatmap <- renderPlot({
    df <- reactiveData()
    nums <- df[sapply(df, is.numeric)]
    if (ncol(nums) > 1) {
      cor_mat <- round(cor(nums, use = "pairwise.complete.obs"), 2)
      cor_df <- melt(cor_mat)
      ggplot(cor_df, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        geom_text(aes(label = value), color = "white") +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
        theme_minimal() +
        labs(title = "Correlation Heatmap", x = "", y = "")
    } else {
      plot.new()
      text(0.5, 0.5, "Not enough numeric variables for correlation heatmap.")
    }
  })
  
  ## Simulation Plot (Cox PH Model)
  output$simulationPlot <- renderPlot({
    df <- reactiveData()
    fit_cox <- coxph(modelFormula(), data = df, x = TRUE, y = TRUE)
    newdata_sim <- getDifferentiatedData(df, input$covariates)
    newdata_sim$treatment <- factor(c("CABG", "PCI"), levels = levels(df$treatment))
    if ("age" %in% input$covariates) newdata_sim$age <- c(input$simAge, input$simAge)
    surv_sim <- tryCatch(survfit(fit_cox, newdata = newdata_sim),
                         error = function(e) { return(NULL) })
    if (is.null(surv_sim) || (is.null(surv_sim$strata) && length(surv_sim$time) == 0)) {
      plot.new()
      text(0.5, 0.5, "Simulation analysis not available for this model.")
    } else {
      ggsurv <- tryCatch({
        ggsurvplot(surv_sim, data = df, risk.table = TRUE, pval = TRUE,
                   legend.title = "Treatment", legend.labs = c("CABG", "PCI"),
                   title = paste("Simulated Survival Curves (Age =", input$simAge, ")"),
                   ggtheme = theme_minimal())
      }, error = function(e) { return(NULL) })
      if (is.null(ggsurv)) {
        plot.new()
        text(0.5, 0.5, "Simulation plot failed.")
      } else {
        print(ggsurv$plot)
      }
    }
  })
  
  ## Validation: Bootstrap for Cox-based models
  output$validationResults <- renderPrint({
    df <- reactiveData()
    modType <- input$modelType
    if (modType %in% c("Cox Proportional Hazards", "Cox with Splines")) {
      dd <- datadist(df)
      options(datadist = dd)
      fit_rms <- cph(modelFormula(), data = df, x = TRUE, y = TRUE, surv = TRUE)
      val <- validate(fit_rms, method = "boot", B = input$nBootstrap)
      print(val)
    } else {
      cat("Bootstrap validation is implemented only for Cox-based models.\n")
    }
  })
  
  ## Extra UI Outputs (Age Histogram, Time Boxplot, Dataset Summary)
  output$extraUI <- renderUI({
    if (input$extraOutputType == "Age Histogram") {
      plotOutput("ageHist")
    } else if (input$extraOutputType == "Time Boxplot") {
      plotOutput("timeBoxplot")
    } else if (input$extraOutputType == "Dataset Summary") {
      verbatimTextOutput("dataSummary")
    }
  })
  
  output$ageHist <- renderPlot({
    df <- reactiveData()
    hist(df$age, main = "Age Histogram", xlab = "Age", col = "lightblue")
  })
  
  output$timeBoxplot <- renderPlot({
    df <- reactiveData()
    boxplot(df$time, main = "Time Boxplot", ylab = "Time")
  })
  
  output$dataSummary <- renderPrint({
    df <- reactiveData()
    summary(df)
  })
}

shinyApp(ui = ui, server = server)
