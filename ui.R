library(shiny)
library(shinyIncubator)
library(rgl)

# Define UI for application that draws a histogram
shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("Understanding Structural Coexistence for Multispecies Communities"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        numericInput("r1","r1",1,min=0,max=1,step=0.1),
        numericInput("r2","r2",1,min=0,max=1,step=0.1),
        numericInput("r3","r3",0.5,min=0,max=1,step=0.1),
        matrixInput("alphamat",
                    "Competition coefficients matrix (α)",
                    data.frame(matrix(c(1,0.5,0.05,0.4,1,0.5,0.3,0.6,1),nrow=3))),
        checkboxInput("draw_table","Show table",value=T),
        radioButtons("cone_opt", "Plot cone",
                     choices = c("none", "static", "interactive"),
                     selected = "none")
        # checkboxInput("draw_cone","Plot cone",value=F),
        # checkboxInput("draw_3Dcone","Plot interactive cone",value=F)
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
        conditionalPanel("input.draw_table",fluidRow(tableOutput("stats")),align="center"),
        plotOutput("proj"),
        conditionalPanel("input.cone_opt == 'static'",
                         plotOutput("cone"), align = "center"),
        conditionalPanel("input.cone_opt == 'interactive'",
                         rglwidgetOutput("cone3d"), align = "center")
      )
    )
  ))

