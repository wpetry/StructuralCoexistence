library(shiny)
library(shinyMatrix)
library(rgl)
library(scales)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel(
      h1("How is diversity maintained in multispecies communities?", h2("A structural approach to coexistence"))
      ),
    
    # tags
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    
    # Inputs
    
    sidebarLayout(
     
      sidebarPanel(
        p("Coexistence theory has largely focused on 2-species communities, but this neglects the role of indirect interactions that are found only in larger systems. Saavedra et al. (doi:10.1002/ecm.1263) offers an approach to quantifying diversity in these systems using the mathematics of structural stability borrowed from engineering. Use the controls below to explore how competitor fitness & interactions govern diversity."),
        radioButtons("spp", "Number of Species",
                     choices = c('3','4'),
                     selected = '3'),
        conditionalPanel(
          condition = "input.spp == '3'",
          h4("Species intrinsic growth rates:"),
          numericInput("r1","r1",1,min=0,max=1,step=0.1),
          numericInput("r2","r2",0.6,min=0,max=1,step=0.1),
          numericInput("r3","r3",0.7,min=0,max=1,step=0.1),
          h4("Interaction coefficients:"),
          matrixInput("alphamat",
                      matrix(c(1, 0.4, 0.3,
                               0.5, 1, 0.6,
                               0.05, 0.5, 1),
                             nrow = 3,
                             dimnames = list(c("α(1,_)", "α(2,_)", "α(3,_)"), c("α(_,1)", "α(_,2)", "α(_,3)")),
                             byrow = TRUE),
                      class = "numeric",
                      rows = list(names = TRUE,
                                  editableNames = FALSE),
                      cols = list(names = TRUE,
                                  editableNames = FALSE)),
          h4("Answers"),
          actionButton("neutral", "Q1. Quasi-neutrality"),
          actionButton("intransient", "Q2. Rock-paper-scissors"),
          actionButton("weakintra", "Q3. Weak interactions"),
          hr(),
          h4("Display options"),
          checkboxInput("draw_table","Show table",value=T),
          checkboxInput("draw_network","Show network",value=T),
          radioButtons("cone_opt", "Plot cone",
                       choices = c("none", "static",
                                   "interactive"),
                       selected = "none")
        ),
        conditionalPanel(
          condition = "input.spp == '4'",
          numericInput("rr1","r1",1,min=0,max=1,step=0.1),
          numericInput("rr2","r2",1,min=0,max=1,step=0.1),
          numericInput("rr3","r3",0.5,min=0,max=1,step=0.1),
          numericInput("rr4","r4",0.5,min=0,max=1,step=0.1),
          matrixInput("alphamat4",
                      #"Competition coefficients matrix (α)",
                      matrix(c(1, 0.4, 0.3, 0.1,
                               0.5, 1, 0.6, 0.1,
                               0.05, 0.5, 1, 0.5,
                               0.1, 0.2, 0.2, 1),
                             byrow = TRUE, nrow = 4)),
          checkboxInput("draw_table4sp","Show table",value=T)
          
        )
   
       ),
    
      # Outputs
      
     mainPanel(
       conditionalPanel("input.draw_table & input.draw_network",
                        fluidRow(
                          column(width = 8,
                                 plotOutput("network")),
                          column(width = 4,
                                 h5("Coexistence metrics"),
                                 tableOutput("stats"))
                        ),
                        align = "center"),
       conditionalPanel("input.spp == '3'",
                        plotOutput("proj"),
                        align = "right"),
       conditionalPanel("input.spp == '3' & input.cone_opt == 'static'",
                        plotOutput("cone"), align = "center"),
       conditionalPanel("input.spp == '3' & input.cone_opt == 'interactive'",
                        rglwidgetOutput("cone3d"), align = "center"),
       conditionalPanel("input.spp == '4' & input.draw_table4sp",
                        fluidRow(tableOutput("stats4sp")),
                        align = "center"),
       conditionalPanel("input.spp == '4'",
                        rglwidgetOutput("proj4sp"), align = "center")
     )
     
    )
  ))

