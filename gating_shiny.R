library(ggplot2)
library(dplyr)
library(shiny)
library(shinyjs)
library(DT)

source('gating_functions.R')

ui <- fluidPage(
  shinyjs::useShinyjs(),
  tags$head(
    tags$style(HTML("
      .custom-title {
        font-size: 36px; 
        font-weight: bold;
        text-align: center; 
        margin-bottom: 20px; 
      }
      .dataTables_wrapper .dataTable {
        width: 100%;
        margin: auto;
        border: 2px solid black;  
      }
      table.dataTable {
        border-collapse: separate !important;
        border-spacing: 0 !important;
      }
      table.dataTable th, table.dataTable td {
        border-left: 1px solid black !important;
        border-top: 1px solid black !important;
      }
      table.dataTable th:first-child, table.dataTable td:first-child {
        border-left: none !important;
      }
      table.dataTable th, table.dataTable td, table.dataTable tr:last-child {
        border-bottom: 1px solid black !important;
      }
      table.dataTable th:last-child, table.dataTable td:last-child {
        border-right: 1px solid black !important;
      }
       .styled-text {
        font-size: 16px; 
        border: 2px solid blue; 
        background-color: #f0f8ff; 
        padding: 10px;
       }
    }
    "))
  ),
  titlePanel(
    tags$div("Calculation of Gating Probabilities", class = "custom-title")
  ),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      hr(style = "border-top: 1px solid #000000;margin-top:0px;margin-bottom:0px;"),
      h2(style = "font-size: 18px; font-weight: bold;margin-top:4px;margin-bottom:4px;", 'Target ORRs'),
      sliderInput("orr_minimum_go", "Minimum ORR Rate (%):", min = 0, max = 1, value = 0.3),
      hr(style = "border-top: 1px solid #000000;margin-top:2px;margin-bottom:2px;"),
      h2(style = "font-size: 18px; font-weight: bold;margin-top:4px;margin-bottom:4px;", 'Probability Thresholds'),
      sliderInput("confidence_go", "GO Confidence (%):", min = 0, max = 1, value = 0.8),
      div(style = "margin-top:-14px"),
      sliderInput("confidence_stop", "STOP Confidence (%):", min = 0, max = 1, value = 0.2),
      hr(style = "border-top: 1px solid #000000;margin-top:2px;margin-bottom:2px;"),
      h2(style = "font-size: 18px; font-weight: bold;margin-top:4px;margin-bottom:4px;", 'Probabilities of Interest'),
      sliderInput("true_orrs", "True ORRs (%):", min = 0, max = 1, value = c(0.1,0.6)),
      hr(style = "border-top: 1px solid #000000;margin-top:2px;margin-bottom:2px;"),
      h2(style = "font-size: 18px; font-weight: bold;margin-top:4px;margin-bottom:-8px;", 'Sample Size'),
      radioButtons("radioGroup", "", 
                   choices = list("Number of Patients (single value)" = "number_pats", 
                                  "Number of Patients (range)" = "number_pats_range"),
                   selected = "number_pats"),
      div(id = "slider1_div", 
          style = "margin-top:-4px;margin-bottom:-6px;", 
          sliderInput("number_pats", "Number of Patients (single value)", min = 0, max = 100, value = 40)),
      div(id = "slider2_div", 
          style = "margin-top:-4px;margin-bottom:-6px;", 
          sliderInput("number_pats_range", "Number of Patients (range)", min = 0, max = 100, value = c(20,60))),
      hr(style = "border-top: 1px solid #000000;margin-top:2px;margin-bottom:2px;")
      
    ),
    mainPanel(
      shinyjs::hidden(
        div(id = "main_content",
            h3("Gating Rules"),
            tags$div(textOutput("go_decision_rule"), style = "font-size: 17px;"),
            tags$div(textOutput("stop_decision_rule"), style = "font-size: 17px;"),
            tags$br(),
            h3("Gating Decision Table"),
            div(id = "multiple_sample_size_section",
                div(
                  #style = "display: flex; align-items: left;",
                  div(style = "margin-top:0px; margin-bottom:-4px; margin-right: 1px; 
                      margin-left: 1px; padding: 0;", 
                      checkboxInput('grouping_option', 'Group cases', value = FALSE)),
                  div(style = "margin-top:0px; margin-bottom:-4px; margin-right: 1px; 
                      margin-left: 1px; padding: 0;", 
                      checkboxInput('stop_option', 'Show stop case', value = TRUE)),
                  div(style = "margin-top:0px; margin-bottom:-4px; margin-right: 1px; 
                      margin-left: 1px; padding: 0;", 
                      checkboxInput('eval_option', 'Show evaluate case', value = TRUE)),
                  div(style = "margin-top:0px; margin-bottom:-4px; margin-right: 1px; 
                      margin-left: 1px; padding: 0;", 
                      checkboxInput('go_option', 'Show go case', value = TRUE))
                )),
            DTOutput("gating_table"),
            tags$br(),
            uiOutput("cell_text"),
            tags$br(),
            div(id = "decision_plot_section",
                h3("Cumulative Decision Probability Plot"),
                checkboxInput('cdpp_show_line', 'Show vertical line', value = TRUE),
                numericInput('cdpp_intercept', 'Choose x-intercept', value = 0.3, 
                             min = 0, max = 1, step = 0.1),
                plotOutput("cumulative_decision_prob_plot"),
                tags$br(),
                h3("Decision Probability Plot"),
                checkboxInput('dpp_show_line', 'Show vertical line', value = TRUE),
                numericInput('dpp_intercept', 'Choose x-intercept', value = 0.3, 
                             min = 0, max = 1, step = 0.1),
                plotOutput("decision_prob_plot")
            )  
        )
      )
    )
  )
)

# Server.
server <- function(input, output, session) {
  
  observeEvent(input$radioGroup, {
    if (input$radioGroup == "number_pats") {
      shinyjs::show("slider1_div")
      shinyjs::hide("slider2_div")
      shinyjs::show("decision_plot_section")
    } else if (input$radioGroup == "number_pats_range") {
      shinyjs::show("slider2_div")
      shinyjs::hide("slider1_div")
      shinyjs::hide("decision_plot_section")
    }
  })
  
  plot_data <- reactive({
    orr_minimum_go <- input$orr_minimum_go
    confidence_go <- input$confidence_go
    confidence_stop <- input$confidence_stop
    prob0 <- input$true_orrs[1]
    prob1 <- input$true_orrs[2]
    number_pats <- input$number_pats
    
    resp_go <- get_responder_go(number_pats, orr_minimum_go, confidence_go)
    validate(
      need(!is_not_valid_number(resp_go), 
           "Invalid number for responders in 'go' decision. Please check your input settings.")
    )
    resp_stop <- get_responder_stop(number_pats, orr_minimum_go, confidence_stop)
    validate(
      need(!is_not_valid_number(resp_stop), 
           "Invalid number for responders in 'stop' decision. Please check your input settings.")
    )
    
    create_oc_plot_data(number_pats, orr_minimum_go, confidence_go, confidence_stop)
  })
  
  cumulative_decision_prob_plot <- reactive({
    
    if (input$radioGroup == "number_pats_range") {
      return(NULL)
    }  
    
    req(plot_data())
    
    intercept <- input$cdpp_intercept
    
    # Calculate intersection points
    go_line <- plot_data() %>% dplyr::filter(group == 1)
    stop_line <- plot_data() %>% dplyr::filter(group == 3)
    evaluate_line <- plot_data() %>% dplyr::filter(group == 2)
    
    go_y <- approx(go_line$prob, go_line$prob.p, xout = intercept)$y
    stop_y <- approx(stop_line$prob, stop_line$prob.p, xout = intercept)$y
    evaluate_y <- approx(evaluate_line$prob, evaluate_line$prob.p, xout = intercept)$y
    
    delta_x <- 0.025 # 0.03
    delta_y <- 0.02 # 0.02
    stop_y_text <- stop_y / 2 + delta_y
    evaluate_y_text <- evaluate_y / 2 + stop_y + delta_y
    go_y_text <- go_y / 2 + evaluate_y + stop_y + delta_y
    
    font_size = 5
    
    plot <- draw_cumulative_decision_probabilities(plot_data())
    
    if (input$cdpp_show_line) {
      plot + 
        geom_vline(xintercept = intercept, linetype = "dashed", color = "black", 
                   linewidth = 0.8) +
        annotate("text", x = intercept - delta_x, y = stop_y_text, 
                 label = paste0('s=',round(stop_y, 2)), color = "black", size = font_size) +
        annotate("text", x = intercept - delta_x, y = evaluate_y_text+0.05, 
                 label = paste0('e=',round(evaluate_y, 2)), color = "black", size = font_size) +
        annotate("text", x = intercept - delta_x, y = go_y_text, 
                 label = paste0('g=',round(go_y, 2)), color = "black", size = font_size) 
    } else {
      draw_cumulative_decision_probabilities(plot_data())
    }
  })
  
  decision_prob_plot <- reactive({
    
    if (input$radioGroup == "number_pats_range") {
      return(NULL)
    }  
    
    req(plot_data())
    
    intercept <- input$dpp_intercept
    
    # Calculate intersection points.
    go_line <- plot_data() %>% dplyr::filter(group == 1)
    stop_line <- plot_data() %>% dplyr::filter(group == 3)
    evaluate_line <- plot_data() %>% dplyr::filter(group == 2)
    
    go_y <- approx(go_line$prob, go_line$prob.p, xout = intercept)$y
    stop_y <- approx(stop_line$prob, stop_line$prob.p, xout = intercept)$y
    evaluate_y <- approx(evaluate_line$prob, evaluate_line$prob.p, xout = intercept)$y
    
    delta_x <- 0.03
    delta_y <- 0.02
    
    font_size = 5
    
    plot <- draw_decision_probabilities(plot_data())
    if (input$dpp_show_line) {
      plot + 
        geom_vline(xintercept = intercept, linetype = "dashed", color = "black", linewidth = 0.8) +
        geom_hline(yintercept = go_y, linetype = "dotted", color = "darkgreen", linewidth = 0.9) +
        geom_hline(yintercept = stop_y, linetype = "dotted", color = "red2", linewidth = 0.9) +
        geom_hline(yintercept = evaluate_y, linetype = "dotted", color = "orange1", linewidth = 0.9) +
        annotate("point", x = intercept, y = go_y, color = "darkgreen", shape = 16, size = 3.5) +
        annotate("text", x = intercept+delta_x, y = go_y+delta_y, 
                 label = paste0("g=", round(go_y, 2)), color = "darkgreen", size = font_size) +
        annotate("point", x = intercept, y = stop_y, color = "red2", shape = 16, size = 3.5) +
        annotate("text", x = intercept+delta_x, y = stop_y-delta_y, 
                 label = paste0("s=", round(stop_y, 2)), color = "red2", size = font_size) +
        annotate("point", x = intercept, y = evaluate_y, color = "orange1", shape = 16, size = 3.5) +
        annotate("text", x = intercept+delta_x, y = evaluate_y+0.05, 
                 label = paste0("e=", round(evaluate_y, 2)), color = "orange1", size = font_size)  
    } else {
      draw_decision_probabilities(plot_data())
    }
  })
  
  go_threshold <- reactive({
    req(plot_data())
    input$orr_minimum_go
  })
  
  stop_threshold <- reactive({
    req(plot_data())
    input$orr_minimum_go
  })
  
  go_probability <- reactive({
    req(plot_data())
    input$confidence_go
  })
  
  stop_probability <- reactive({
    req(plot_data())
    input$confidence_stop
  })
  
  gating_table <- reactive({
    
    req(plot_data())
    
    if (input$radioGroup == "number_pats") {
      output$cell_text <- renderUI({ NULL })
      make_decision_table_for_prob_of_interest(input$number_pats, 
                                               input$orr_minimum_go,
                                               input$confidence_go,
                                               input$confidence_stop,
                                               input$true_orrs,
                                               input$stop_option,
                                               input$go_option,
                                               input$eval_option)   
    } else {
      output$cell_text <- renderUI({ NULL })
      grouping_option <- input$grouping_option
      stop_option <- input$stop_option
      go_option <- input$go_option
      eval_option <- input$eval_option
      make_decision_table_for_multiple_sample_sizes(input$number_pats_range, 
                                                    input$orr_minimum_go,
                                                    input$confidence_go,
                                                    input$confidence_stop,
                                                    input$true_orrs,
                                                    input$grouping_option,
                                                    input$stop_option,
                                                    input$go_option,
                                                    input$eval_option)
    }  
  })
  
  output$cumulative_decision_prob_plot <- renderPlot({
    cumulative_decision_prob_plot()
  })
  
  output$decision_prob_plot <- renderPlot({
    decision_prob_plot()
  })
  
  output$go_decision_rule <- renderText({
    go_threshold <- go_threshold()
    go_threshold <-  as.numeric(go_threshold) * 100
    go_probability <- go_probability()
    go_probability <-  as.numeric(go_probability) * 100
    paste0("Gating rule for GO: P[ORR ≥ ", 
           go_threshold, 
           "% | observed data] ≥ ", 
           go_probability, "%")
  })
  
  output$stop_decision_rule <- renderText({
    stop_threshold <- stop_threshold()
    stop_threshold <-  as.numeric(stop_threshold) * 100
    stop_probability <- stop_probability()
    stop_probability <-  as.numeric(stop_probability) * 100
    paste0("Gating rule for STOP: P[ORR ≥ ",
           stop_threshold, 
           "% | observed data] < ", 
           stop_probability, "%")
  })
  
  output$gating_table <- renderDT({
    
    req(gating_table())
    
    df <- gating_table()
    num_cols <- ncol(df)
    if (input$radioGroup == "number_pats") {
      DT::datatable(df, 
                    colnames = rep("", ncol(df)),
                    rownames = FALSE,
                    selection = 'single', 
                    extensions = 'Buttons',
                    options = list(
                      dom = 'Bt', 
                      buttons = list(
                        list(extend = 'copy', title = "Gating Decision Table Data"),
                        list(extend = 'csv', title = "Gating Decision Table Data"),
                        list(extend = 'excel', title = "Gating Decision Table Data")
                      ),
                      paging = FALSE,
                      columnDefs = list(
                        list(width = '20%', targets = 0),
                        list(width = '30%', targets = 1)
                      ),
                      rowCallback = JS(
                        "function(row, data, displayNum, index){",
                        "if(data[0] == 'True ORR'){",
                        "$(row).css('background-color', '#ADD8E6');",
                        "$(row).css('color', 'white');",
                        "}",
                        "}"
                      ),
                      headerCallback = JS(
                        "function(thead, data, start, end, display){",
                        "  $(thead).remove();",
                        "}"
                      )
                    ),
                    class = "display") %>%
        DT::formatStyle(
          columns = c(1,2,3,4),
          backgroundColor = styleEqual(
            c("STOP", "EVALUATE", "GO"),
            c("#EE0000", "#FFA500", "darkgreen")
          )
        )  
    } else {
      output$cell_text <- renderUI({ NULL })
      DT::datatable(df, 
                    colnames = rep("", num_cols),
                    rownames = FALSE,
                    selection = 'single', 
                    extensions = 'Buttons',
                    options = list(
                      dom = 'Bt', 
                      buttons = list(
                        list(extend = 'copy', title = "Gating Decision Table Data"),
                        list(extend = 'csv', title = "Gating Decision Table Data"),
                        list(extend = 'excel', title = "Gating Decision Table Data")
                      ),
                      paging = FALSE,
                      columnDefs = list(
                        list(width = '20%', targets = 0),
                        list(width = '30%', targets = 1)
                      ),
                      rowCallback = JS(
                        "function(row, data, displayNum, index){",
                        "if(data[0] == 'True ORR'){",
                        "$(row).css('background-color', '#ADD8E6');",
                        "$(row).css('color', 'white');",
                        "}",
                        "}"
                      ),
                      headerCallback = JS(
                        "function(thead, data, start, end, display){",
                        "  $(thead).remove();",
                        "}"
                      )
                    ),
                    class = "display") %>%
        DT::formatStyle(
          columns = seq(2, num_cols),  # Dynamically generate column indices
          backgroundColor = styleEqual(
            c("STOP", "EVALUATE", "GO"),
            c("#EE0000", "#FFA500", "darkgreen")
          )
        )  
    }
  })
  
  # Display text when a row is clicked.
  observeEvent(input$gating_table_rows_selected, {
    
    row <- input$gating_table_rows_selected
    df <- gating_table()
    selected_values <- df[row, , drop = FALSE]
    probs <- unlist(selected_values)
    res <- grepl("%$", probs)
    
    if (ncol(df) != 4) {
      output$cell_text <- renderUI({ NULL })
      return(NULL)
    }
    if (sum(res) == 4 && all(res == TRUE)) {
      true_orr <- probs[1]
      stop_prob <- probs[2]
      eval_prob <- probs[3]
      go_prob <- probs[4]
      obs_orr_stop <- df[3,2]
      obs_orr_eval <- df[3,3]
      obs_orr_go <- df[3,4]
      target_orr_percent <- as.numeric(input$orr_minimum_go) * 100
      true_orr_percent <- as.numeric(gsub("%", "", true_orr))
      
      t1 <- paste0('Assuming a true ORR of ',true_orr,' and based on the observed data with ORR ranges of ',
                   obs_orr_go,' for a GO decision and ',obs_orr_stop,
                   ' for a STOP decision, we have the following probabilities:')
      # True GO: true ORR above target ORR.
      # Wrong GO: true ORR below target ORR.
      # True STOP: true ORR below target ORR.
      # Wrong STOP: true ORR above target ORR.
      if (true_orr_percent >= target_orr_percent) {
        t2 <- paste0('There is a ',
                     go_prob,
                     ' probability of a true GO decision (i.e. correctly proceeding when the true ORR is equal or above the target ORR of ',
                     target_orr_percent,'%)')
      } else {
        t2 <- paste0('There is a ',
                     go_prob,
                     ' probability of a wrong GO decision (i.e. falsely proceeding when the true ORR is below the target ORR of ',
                     target_orr_percent,'%)')
      }
      if (true_orr_percent < target_orr_percent) {
        t3 <- paste0('There is a ',
                     stop_prob,
                     ' probability of a true STOP decision (i.e. correctly stopping when the true ORR is below the target ORR of ',
                     target_orr_percent,'%)')
      } else {
        t3 <- paste0('There is a ',
                     stop_prob,
                     ' probability of a wrong STOP decision (i.e. falsely stopping when the true ORR is equal or above the target ORR of ',
                     target_orr_percent,'%)')
      }  
      t4 <- paste0('There is a ',
                   eval_prob,
                   ' probability of needing further evaluation (i.e. the decision is neither clearly GO nor STOP)')
      
      output$cell_text <- renderUI({
        tags$div(
          class = "styled-text",
          tags$p(t1),
          tags$ul(
            tags$li(t2),
            tags$li(t3),
            tags$li(t4)
          )
        )
      })
      
      output$cell_text <- renderUI({
        tags$div(
          class = "styled-text",
          tags$p(t1),
          tags$ul(
            tags$li(t2),
            tags$li(t3),
            tags$li(t4)
          )
        )
      })
    } else {
      output$cell_text <- renderUI({ NULL })
    }
  })
  
  observeEvent({
    input$orr_minimum_go
    input$confidence_go
    input$confidence_stop
    input$true_orrs
    input$number_pats
  }, {
    output$cell_text <- renderUI({
      NULL
    })
  })
  
  observeEvent(plot_data(), {
    req(gating_table())
    req(cumulative_decision_prob_plot())
    shinyjs::show(id = "main_content")
  }, once = TRUE)
  
}

shinyApp(ui, server)