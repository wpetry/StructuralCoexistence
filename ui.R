#################################################-
## USER INTERFACE
#################################################-
## Preliminaries ----
#################################################-
library(shiny)
library(shinyMatrix)
library(rgl)
library(scales)
#################################################-
## Define UI ----
#################################################-
shinyUI(
  fluidPage(
    # Application title
    titlePanel(
      h1("How is diversity maintained in multispecies communities?", h2("A structural approach to coexistence"))
    ),
    # HTML tags
    tags$head(
      tags$style(HTML("hr {border-top: 1px solid #000000;}"))
    ),
    # User inputs
    sidebarLayout(
      sidebarPanel(
        p("Coexistence theory has largely focused on 2-species communities, but this neglects the role of indirect interactions that are found only in larger systems. Saavedra et al. (doi:10.1002/ecm.1263) offers an approach to quantifying diversity in these systems using the mathematics of structural stability borrowed from engineering. Use the controls below to explore how competitor fitness & interactions govern diversity."),
        hr(),
        h4("Number of Species:"),
        radioButtons("spp", label = NULL,
                     choices = c('3', '4'),
                     selected = '3'),
        hr(),
        h4("Species intrinsic growth rates:"),
        conditionalPanel(
          condition = "input.spp == '3'",
          
          numericInput("r1", "r1", 1, min = 0, max = 1,
                       step = 0.1),
          numericInput("r2", "r2", 0.6, min = 0, max = 1,
                       step = 0.1),
          numericInput("r3", "r3", 0.7, min = 0, max = 1,
                       step = 0.1),
          h4("Interaction coefficients:"),
          shinyMatrix::matrixInput(inputId = "alphamat",
                                   value = matrix(c(1, 0.4, 0.3,
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
          hr(),
          conditionalPanel(condition = "input.spp == '3'",
                           h4("Scenarios:"),
                           actionButton("neutral3",
                                        "Quasi-neutrality"),
                           actionButton("intransient3",
                                        "Rock-paper-scissors"),
                           actionButton("weakintra3",
                                        "Weak inter-specific interactions")
          ),
          hr(),
          h4("Display options:"),
          checkboxInput("draw_network", "Species interaction network", value = TRUE),
          checkboxInput("draw_table", "Coexistence metric table", value = TRUE),
          radioButtons("cone_opt", "Plot cone",
                       choices = c("none", "static", "interactive"),
                       selected = "none")
        ),
        conditionalPanel(
          condition = "input.spp == '4'",
          numericInput("rr1", "r1", 1, min = 0, max = 1, step = 0.1),
          numericInput("rr2", "r2", 1, min = 0, max = 1, step = 0.1),
          numericInput("rr3", "r3", 0.5, min = 0, max = 1, step = 0.1),
          numericInput("rr4", "r4", 0.5, min = 0, max = 1, step = 0.1),
          h4("Interaction coefficients:"),
          shinyMatrix::matrixInput(inputId = "alphamat4",
                                   value = matrix(c(1, 0.4, 0.3, 0.1,
                                                    0.5, 1, 0.6, 0.1,
                                                    0.05, 0.5, 1, 0.5,
                                                    0.1, 0.2, 0.2, 1),
                                                  nrow = 4,
                                                  dimnames = list(c("α(1,_)", "α(2,_)", "α(3,_)", "α(4,_)"),
                                                                  c("α(_,1)", "α(_,2)", "α(_,3)", "α(_,4)")),
                                                  byrow = TRUE),
                                   class = "numeric",
                                   rows = list(names = TRUE,
                                               editableNames = FALSE),
                                   cols = list(names = TRUE,
                                               editableNames = FALSE)),
          hr(),
          conditionalPanel(condition = "input.spp == '4'",
                           h4("Scenarios:"),
                           actionButton("neutral4",
                                        "Quasi-neutrality"),
                           actionButton("weakintra4",
                                        "Weak inter-specific interactions")
          ),
          hr(),
          h4("Display options:"),
          checkboxInput("draw_network4sp", "Species interaction network", value = TRUE),
          checkboxInput("draw_table4sp", "Coexistence metric table", value = TRUE)
        )
      ),
      
      # Outputs
      mainPanel(
        # 3 species
        fluidRow(
          conditionalPanel("input.spp == '3' & input.draw_network",
                           column(width = 7,
                                  plotOutput("network_3sp")
                           )
          ),
          conditionalPanel("input.spp == '3' & input.draw_table",
                           column(width = 5,
                                  h5("Coexistence metrics"),
                                  tableOutput("stats")
                           )
          )),
        conditionalPanel("input.spp == '3'",
                         fluidRow(
                           plotOutput("proj")
                         )
        ),
        conditionalPanel("input.spp == '3' & input.cone_opt == 'static'",
                         fluidRow(
                           plotOutput("cone")
                         )
        ),
        conditionalPanel("input.spp == '3' & input.cone_opt == 'interactive'",
                         fluidRow(
                           rglwidgetOutput("cone3d")
                         )
        ),
        # 4 species
        fluidRow(
          conditionalPanel("input.spp == '4' & input.draw_network4sp",
                           column(width = 7,
                                  plotOutput("network_4sp")
                           )
          ),
          conditionalPanel("input.spp == '4' & input.draw_table4sp",
                           column(width = 5,
                                  h5("Coexistence metrics"),
                                  tableOutput("stats4sp")
                           )
          )),
        conditionalPanel("input.spp == '4'",
                         fluidRow(
                           rglwidgetOutput("proj4sp")
                         )
        )
      )
    )
  )
)
